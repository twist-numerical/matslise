import Problem from "./Problem";

export default class MatsliseController {
  public eigenvalues: [number, number, number][] | null = null;
  public problem = new Problem();

  reset() {
    this.problem.reset();
    this.eigenvalues = null;
  }

  calculate() {
    this.reset();
    this.eigenvalues = [];
    this.moreEigenvalues();
  }

  moreEigenvalues() {
    for (const [index, value] of this.problem.eigenvaluesByIndex(
      this.eigenvalues!.length,
      this.eigenvalues!.length + 10
    )) {
      this.eigenvalues!.push([
        index,
        value,
        this.problem.eigenvalueError(index, value)
      ]);
    }
  }
}
