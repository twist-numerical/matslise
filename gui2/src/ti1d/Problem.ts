import Module, {
  AbstractMatslise,
  MatsliseHalfRange,
  Matslise
} from "../util/matsliseLoader";
import math from "mathjs-expression-parser";

const evaluatePair = ([a, b]: [string, string]): [number, number] => {
  return [math.eval(a), math.eval(b)];
};

export default class Problem {
  public Matslise: typeof Matslise | null = null;
  public HalfRange: typeof MatsliseHalfRange | null = null;
  public matslise: AbstractMatslise | null = null;
  public parsed: {
    potential: (x: number) => number;
    x: [number, number];
    ymin: [number, number];
    ymax: [number, number];
    left: [number, number];
    right: [number, number];
    tolerance: number;
    symmetric: boolean;
  } | null = null;

  public potential = "0";
  public x: [string, string] = ["0", "pi"];
  public ymin: [string, string] = ["1", "0"];
  public ymax: [string, string] = ["1", "0"];
  public tolerance = "1e-8";
  public symmetric = false;
  public toDelete: { delete: () => void }[] = [];

  constructor() {
    this.initMatslise();
  }

  initMatslise() {
    new Module().then(module => {
      this.Matslise = module.Matslise;
      this.HalfRange = module.MatsliseHalfRange;
    });
  }

  parse() {
    if (this.Matslise === null || this.HalfRange === null)
      throw new Error("Wait until problem.Matslise !== null");
    if (this.matslise !== undefined) this.reset();
    const compiledPotential = math.compile(this.potential);
    const potential = (x: number) => compiledPotential.eval({ x });

    const ymin = evaluatePair(this.ymin);
    const ymax = evaluatePair(this.ymax);
    const x = evaluatePair(this.x);
    this.parsed = {
      potential,
      x: this.symmetric ? [-x[1], x[1]] : x,
      ymin: this.symmetric ? ymax : ymin,
      ymax,
      left: this.symmetric ? [ymax[1], -ymax[0]] : [ymin[1], -ymin[0]],
      right: [ymax[1], -ymax[0]],
      tolerance: math.eval(this.tolerance),
      symmetric: this.symmetric
    };
    this.matslise = this.parsed.symmetric
      ? new this.HalfRange(potential, this.parsed.x[1], {
          tolerance: this.parsed.tolerance
        })
      : new this.Matslise(potential, this.parsed.x[0], this.parsed.x[1], {
          tolerance: this.parsed.tolerance
        });
    this.toDelete.push(this.matslise!);
  }

  reset() {
    this.parsed = null;
    this.matslise = null;
    for (const obj of this.toDelete.splice(0, this.toDelete.length))
      obj.delete();
  }

  eigenvaluesByIndex(imin: number, imax: number): [number, number][] {
    if (this.matslise === null || this.parsed === null) this.parse();
    return this.matslise!.eigenvaluesByIndex(
      imin,
      imax,
      this.parsed!.left,
      this.parsed!.right
    ).map(({ first, second }) => [first, second]);
  }

  eigenvalueError(index: number, E: number): number {
    if (this.matslise === null || this.parsed === null) this.parse();
    return this.matslise!.eigenvalueError(
      E,
      this.parsed!.left,
      this.parsed!.right,
      index
    );
  }

  eigenfunction(index: number, E: number, x: number[]): number[] {
    if (this.matslise === null || this.parsed === null) this.parse();

    return this.matslise!.computeEigenfunction(
      E,
      this.parsed!.left,
      this.parsed!.right,
      x,
      index
    );
  }
}
