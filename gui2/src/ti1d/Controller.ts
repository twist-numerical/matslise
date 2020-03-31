import Problem from "./Problem";
import Color = require("color");

export interface Eigenvalue {
  index: number;
  value: number;
  error: number;
  visible: boolean;
  color: string;
  readonly eigenfunction: number[];
}

export default class Controller {
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
    const n = 120;
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
    this.eigenvalues[0].visible = true;
    this.eigenvalues[1].visible = true;
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
        color: Color.hsv(1.61803398875 * index * 100, 100, 80).string(),
        get eigenfunction(): number[] {
          if (self.eigenfunctions.has(index))
            return self.eigenfunctions.get(index)!;
          const f = self.problem.eigenfunction(index, value, self.xValues);
          self.eigenfunctions.set(index, f);
          return f;
        }
      });
    }
  }
}
