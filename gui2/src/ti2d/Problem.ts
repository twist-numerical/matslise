// @ts-ignore
import WorkerPromise from "webworker-promise";

interface HistoryEntry {
  name: string;
  time: number | null;
}

export default class Problem {
  public potential = "(1+x^2)*(1+y^2)";
  public x: [string, string] = ["-5.5", "5.5"];
  public y: [string, string] = ["-5.5", "5.5"];
  public tolerance = "1e-5";
  public xSymmetric = true;
  public ySymmetric = false;
  public worker: WorkerPromise;
  public potentialData: number[][];

  public ready: boolean = false;
  public calculating: boolean = false;
  public parsed: null | {
    tolerance: number;
    x: [number, number];
    y: [number, number];
    xPoints: number[];
    yPoints: number[];
    xSymmetric: boolean;
    ySymmetric: boolean;
  } = null;
  public history: HistoryEntry[] = [];

  constructor() {
    this.startWorker();
  }

  startWorker() {
    if (this.worker === undefined) {
      this.worker = new WorkerPromise(new Worker("./worker.ts"));
      this.worker.postMessage({ type: "isReady" }).then((ready: boolean) => {
        this.ready = ready;
      });
    } else {
      console.error("Worker already started");
    }
  }

  terminate() {
    this.worker.terminate();
    this.worker = undefined;
    this.ready = false;
    this.calculating = false;
    this.parsed = null;
    this.startWorker();
  }

  async parse() {
    this.parsed = await this.addToHistory(
      `Initialising problem`,
      this.worker.postMessage({
        type: "parse",
        data: {
          potential: this.potential,
          x: this.x,
          y: this.y,
          tolerance: this.tolerance,
          xSymmetric: this.xSymmetric,
          ySymmetric: this.ySymmetric,
        },
      })
    );
    this.potentialData = await this.worker.postMessage({
      type: "evaluatePotential",
    });
  }

  reset() {
    if (this.calculating) this.terminate();
    this.parsed = null;
  }

  async addToHistory(name: string, action: Promise<any>): Promise<any> {
    const entry: HistoryEntry = {
      name,
      time: null,
    };
    this.history.push(entry);
    const start = +new Date();
    const result = await action;
    entry.time = +new Date() - start;
    return result;
  }

  async eigenvalues(emin: number, emax: number): Promise<[number, number][]> {
    return await this.addToHistory(
      `Eigenvalues in [${emin.toPrecision(3)}, ${emax.toPrecision(3)}]`,
      this.worker.postMessage({
        type: "eigenvalues",
        data: { emin, emax },
      })
    );
  }

  async eigenvaluesByIndex(
    imin: number,
    imax: number
  ): Promise<[number, number][]> {
    return await this.addToHistory(
      `Eigenvalues`,
      this.worker.postMessage({
        type: "eigenvaluesByIndex",
        data: { imin, imax },
      })
    );
  }

  async firstEigenvalue(): Promise<[number, number]> {
    return await this.addToHistory(
      `First eigenvalues`,
      await this.worker.postMessage({
        type: "firstEigenvalue",
      })
    );
  }

  async eigenfunction(E: number): Promise<number[][][]> {
    return await this.addToHistory(
      `Eigenfunction of ${E.toPrecision(4)}`,
      await this.worker.postMessage({
        type: "eigenfunction",
        data: { E },
      })
    );
  }
}
