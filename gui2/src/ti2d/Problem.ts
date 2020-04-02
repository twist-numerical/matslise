import WorkerPromise from "webworker-promise";

interface HistoryEntry {
  name: string;
  time: number | null;
}

export default class Problem {
  public potential = "(1+x*x)*(1+y*y)";
  public x: [string, string] = ["-3", "3"];
  public y: [string, string] = ["-3", "3"];
  public tolerance = "1e-5";
  public xSymmetric = false;
  public ySymmetric = false;
  public worker: WorkerPromise;

  public ready: boolean = false;
  public calculating: boolean = false;
  public parsed: boolean = false;
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
    this.parsed = false;
    this.startWorker();
  }

  async parse() {
    await this.addToHistory(
      `Initialising problem`,
      this.worker.postMessage({
        type: "parse",
        data: {
          potential: this.potential,
          x: this.x,
          y: this.y,
          tolerance: this.tolerance,
          xSymmetric: this.xSymmetric,
          ySymmetric: this.ySymmetric
        }
      })
    );
    this.parsed = true;
  }

  reset() {
    if (this.calculating) this.terminate();
    this.parsed = false;
  }

  async addToHistory(name: string, action: Promise<any>): Promise<any> {
    const entry: HistoryEntry = {
      name,
      time: null
    };
    this.history.push(entry);
    const start = +new Date();
    const result = await action;
    entry.time = +new Date() - start;
    return result;
  }

  async eigenvalues(emin: number, emax: number): Promise<number[]> {
    return await this.addToHistory(
      `Eigenvalues in [${emin.toPrecision(3)}, ${emax.toPrecision(3)}]`,
      this.worker.postMessage({
        type: "eigenvalues",
        data: { emin, emax }
      })
    );
  }

  async eigenvaluesByIndex(imin: number, imax: number): Promise<number[]> {
    return await this.addToHistory(
      `Eigenvalues`,
      this.worker.postMessage({
        type: "eigenvaluesByIndex",
        data: { imin, imax }
      })
    );
  }

  async firstEigenvalue(): Promise<number> {
    return await this.worker.postMessage({
      type: "firstEigenvalue"
    });
  }

  eigenvalueError(E: number): number {}
}
