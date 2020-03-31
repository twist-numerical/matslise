import Problem from "./Problem";
import Color = require("color");

export interface Eigenvalue {
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

  async calculate() {
    this.reset();
    await this.problem.parse();
    await this.moreEigenvalues();
  }

  async moreEigenvalues() {
    let eigenvaluesFound;
    if (this.eigenvalues === null) {
      eigenvaluesFound = [await this.problem.firstEigenvalue()];
      this.eigenvalues = [];
    } else {
      eigenvaluesFound = await this.problem.eigenvalues(0, 10);
    }
    const self = this;
    for (const value of eigenvaluesFound) {
      const index = 0;
      this.eigenvalues!.push({
        value,
        error: this.problem.eigenvalueError(value),
        visible: false,
        color: Color.hsv(1.61803398875 * index * 100, 100, 80).string(),
        get eigenfunction(): number[] {
          if (self.eigenfunctions.has(index))
            return self.eigenfunctions.get(index)!;
          const f: number[] = []; // self.problem.eigenfunction(index, value, self.xValues);
          //self.eigenfunctions.set(index, f);
          return f;
        }
      });
    }
  }
}
