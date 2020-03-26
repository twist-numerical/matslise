declare class AbstractMatslise {
  delete(): void;

  eigenvaluesByIndex(
    imin: number,
    imax: number,
    left: [number, number],
    right: [number, number]
  ): { first: number; second: number }[];

  eigenvalueError(
    E: number,
    left: [number, number],
    right: [number, number],
    index?: number
  ): number;

  computeEigenfunction(
    E: number,
    left: [number, number],
    right: [number, number],
    x: number[],
    index?: number
  ): number[];

  delete(): void;
}

declare class Matslise extends AbstractMatslise {
  constructor(
    potential: (x: number) => number,
    xmin: number,
    xmax: number,
    options:
      | {
          sectorCount: number;
        }
      | {
          tolerance: number;
        }
  );
}

declare class MatsliseHalfRange extends AbstractMatslise {
  constructor(
    potential: (x: number) => number,
    xmax: number,
    options:
      | {
          sectorCount: number;
        }
      | {
          tolerance: number;
        }
  );
}

declare interface matslise {
  Matslise: typeof Matslise;
  MatsliseHalfRange: typeof MatsliseHalfRange;
}

export default class matsliseLoader {
  /*
  class HalfRange {
    eigenvaluesByIndex(
      imin: number,
      imax: number,
      right: [number, number]
    ): [number, number][];
  }*/

  then(resolve: (matslise: matslise) => void): void;
}
