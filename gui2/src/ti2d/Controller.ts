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

  async calculate() {
    this.reset();
    await this.problem.parse();
    await this.moreEigenvalues();
  }

  async moreEigenvalues() {
    let eigenvaluesFound;
    if (this.eigenvalues === null) {
      eigenvaluesFound = await this.problem.eigenvaluesByIndex(0, 2);
      this.eigenvalues = [];
    } else {
      const max = this.eigenvalues.length - 1;
      const min = Math.max(0, max - 5);
      eigenvaluesFound = await this.problem.eigenvalues(
        this.eigenvalues[max].value,
        2 * this.eigenvalues[max].value - this.eigenvalues[min].value
      );
    }
    const self = this;
    for (const value of eigenvaluesFound) {
      const index = this.eigenvalues!.length;
      this.eigenvalues!.push({
        index,
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
