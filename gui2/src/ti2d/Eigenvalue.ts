import Problem from "./Problem";

export default class Eigenvalue {
  public index: number;
  public value: number;
  public error: number;
  public problem: Problem;
  public visible: boolean = false;
  public eigenfunctions: number[][][] | null = null;
  public calculatingEigenfunction: boolean = false;

  constructor(index: number, value: number, error: number, problem: Problem) {
    this.index = index;
    this.value = value;
    this.error = error;
    this.problem = problem;
  }

  public async ensureEigenfunctions() {
    if (!this.calculatingEigenfunction && this.eigenfunctions === null) {
      this.calculatingEigenfunction = true;
      this.eigenfunctions = await this.problem.eigenfunction(this.value);
    }
  }
}
