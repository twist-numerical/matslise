import WorkerPromise from "webworker-promise";

export default class Problem {
  public potential = "(1+x*x)*(1+y*y)";
  public x: [string, string] = ["-3", "3"];
  public y: [string, string] = ["-3", "3"];
  public tolerance = "1e-7";
  public xSymmetric = false;
  public ySymmetric = false;
  public worker: WorkerPromise;
  public ready: boolean = false;

  constructor() {
    this.worker = new WorkerPromise(new Worker("./worker.ts"));
    this.worker.postMessage({ type: "isReady" }).then((ready: boolean) => {
      this.ready = ready;
    });
  }

  async parse() {
    await this.worker.postMessage({
      type: "parse",
      data: {
        potential: this.potential,
        x: this.x,
        y: this.y,
        tolerance: this.tolerance,
        xSymmetric: this.xSymmetric,
        ySymmetric: this.ySymmetric
      }
    });
  }

  reset() {}

  async eigenvalues(emin: number, emax: number): Promise<number[]> {
    return await this.worker.postMessage({
      type: "eigenvalues",
      data: { emin, emax }
    });
  }

  async firstEigenvalue(): Promise<number> {
    return await this.worker.postMessage({
      type: "firstEigenvalue"
    });
  }

  eigenvalueError(E: number): number {}
}
