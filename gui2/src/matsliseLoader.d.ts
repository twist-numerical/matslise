declare class Matslise {
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
    right: [number, number]
  ): number;
}

declare class HalfRange {
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

  delete(): void;

  eigenvaluesByIndex(
    imin: number,
    imax: number,
    right: [number, number]
  ): { first: number; second: number }[];

  eigenvalueError(E: number, right: [number, number], index: number): number;
}

declare interface matslise {
  Matslise: typeof Matslise;
  HalfRange: typeof HalfRange;
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
