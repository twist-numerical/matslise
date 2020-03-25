import Problem from "./Problem";

export interface Eigenvalue {
  index: number;
  value: number;
  error: number;
  visible: boolean;
  readonly eigenfunction: number[];
}

export default class MatsliseController {
  public eigenvalues: Eigenvalue[] | null = null;
  public problem = new Problem();
  public xValues: number[] = [];
  private eigenfunctions: Map<number, number[]> = new Map();

  reset() {
    this.problem.reset();
    this.eigenfunctions.clear();
    this.eigenvalues = null;
  }

  calculate() {
    this.reset();
    this.eigenvalues = [];
    this.problem.parse();
    this.xValues = [];
    const n = 60;
    const [xmin, xmax] = this.problem.parsed!.x;
    for (let i = 0; i < n; ++i) {
      this.xValues.push(
        i === 0
          ? xmin
          : i === n - 1
          ? xmax
          : xmin + ((xmax - xmin) * (i + Math.random())) / n
      );
    }
    this.moreEigenvalues();
    this.eigenvalues[3].visible = true;
  }

  moreEigenvalues() {
    const self = this;
    for (const [index, value] of this.problem.eigenvaluesByIndex(
      this.eigenvalues!.length,
      this.eigenvalues!.length + 10
    )) {
      this.eigenvalues!.push({
        index,
        value,
        error: this.problem.eigenvalueError(index, value),
        visible: false,
        get eigenfunction() {
          if (self.eigenfunctions.has(index))
            return self.eigenfunctions.get(index);
          const f = self.problem.eigenfunction(index, value, self.xValues);
          self.eigenfunctions.set(index, f);
          return f;
        }
      });
    }
  }
}
