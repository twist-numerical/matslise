// @ts-ignore
import registerWebworker from "webworker-promise/lib/register";
import MatsliseModule, { Matslise2D } from "../util/matsliseLoader";
import math from "mathjs-expression-parser";

const evaluatePair = ([a, b]: [string, string]): [number, number] => {
  return [math.eval(a), math.eval(b)];
};

let _Matslise2D: typeof Matslise2D;
const waitForMatslise: ((success: boolean) => void)[] = [];
setTimeout(() => {
  new MatsliseModule({
    locateFile: (path: string) => {
      return path;
    },
  }).then((matslise) => {
    _Matslise2D = matslise.Matslise2D;
    for (const res of waitForMatslise) {
      res(true);
    }
  });
}, 2000);

let parsed:
  | {
      potential: (x: number, y: number) => number;
      x: [number, number];
      y: [number, number];
      xPoints: number[];
      yPoints: number[];
      tolerance: number;
      xSymmetric: boolean;
      ySymmetric: boolean;
    }
  | undefined;
let matslise: Matslise2D | undefined;
const toDelete: { delete(): void }[] = [];

function reset() {
  matslise = undefined;
  parsed = undefined;
  for (let obj of toDelete) obj.delete();
  toDelete.splice(0, toDelete.length);
}

function parse(data: {
  potential: string;
  x: [string, string];
  y: [string, string];
  tolerance: string;
  xSymmetric: boolean;
  ySymmetric: boolean;
}) {
  if (_Matslise2D === null) throw new Error("Wait until matslise is loaded");
  if (matslise !== undefined) reset();
  const compiledPotential = math.compile(data.potential);
  const potential = (x: number, y: number) => compiledPotential.eval({ x, y });

  const x = evaluatePair(data.x);
  const y = evaluatePair(data.y);
  parsed = {
    potential,
    x: data.xSymmetric ? [-x[1], x[1]] : x,
    y: data.ySymmetric ? [-y[1], y[1]] : y,
    tolerance: math.eval(data.tolerance),
    xSymmetric: data.xSymmetric,
    ySymmetric: data.ySymmetric,
    xPoints: [],
    yPoints: [],
  };

  const n = 40;
  const step = Math.max(
    (parsed.x[1] - parsed.x[0]) / n,
    (parsed.y[1] - parsed.y[0]) / n
  );
  for (let x = parsed.x[0] + step / 2; x < parsed.x[1]; x += step)
    parsed.xPoints.push(x);
  for (let y = parsed.y[0] + step / 2; y < parsed.y[1]; y += step)
    parsed.yPoints.push(y);

  matslise = new _Matslise2D(
    potential,
    parsed.x[0],
    parsed.x[1],
    parsed.y[0],
    parsed.y[1],
    {
      tolerance: parsed.tolerance,
      nested: {
        tolerance: parsed.tolerance,
        symmetric: data.xSymmetric,
      },
    }
  );

  const { potential: _, ...parsedValues } = parsed;
  return parsedValues;
}

function addEigenvalueErrors(eigs: number[]): [number, number][] {
  return eigs.map((E) => [E, matslise!.eigenvalueError(E)]);
}

function eigenvalues(emin: number, emax: number): [number, number][] {
  if (matslise === undefined) throw new Error("Problem not parsed");
  return addEigenvalueErrors(matslise.eigenvalues(emin, emax));
}

function eigenvaluesByIndex(imin: number, imax: number): [number, number][] {
  if (matslise === undefined) throw new Error("Problem not parsed");
  return addEigenvalueErrors(matslise.eigenvaluesByIndex(imin, imax));
}

function firstEigenvalue(): [number, number] {
  if (matslise === undefined) throw new Error("Problem not parsed");
  return addEigenvalueErrors([matslise.firstEigenvalue()])[0];
}

function evaluatePotential(): number[][] {
  if (parsed === undefined) throw new Error("Problem not parsed");
  const potential: number[][] = [];
  for (const y of parsed.yPoints) {
    const row: number[] = [];
    for (const x of parsed.xPoints) {
      row.push(parsed.potential(x, y));
    }
    potential.push(row);
  }
  return potential;
}

function eigenfunction(E: number): number[][][] {
  if (matslise === undefined || parsed === undefined)
    throw new Error("Problem not parsed");
  const result = matslise.computeEigenfunction(
    E,
    parsed.xPoints,
    parsed.yPoints
  );
  console.log(result);
  return result;
}

registerWebworker(async (message: { type: string; data: any }) => {
  console.log(message);
  switch (message.type) {
    case "isReady":
      if (_Matslise2D !== undefined) return true;
      return await new Promise((res, rej) => {
        waitForMatslise.push(res);
      });
    case "parse":
      return parse(message.data);
    case "eigenvalues":
      return eigenvalues(message.data.emin, message.data.emax);
    case "eigenvaluesByIndex":
      return eigenvaluesByIndex(message.data.imin, message.data.imax);
    case "firstEigenvalue":
      return firstEigenvalue();
    case "evaluatePotential":
      return evaluatePotential();
    case "eigenfunction":
      return eigenfunction(message.data.E);
  }
});
