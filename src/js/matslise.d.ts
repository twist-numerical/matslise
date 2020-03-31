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

declare class Matslise2D {
    constructor(
        potential: (x: number, y: number) => number,
        xmin: number,
        xmax: number,
        ymin: number,
        ymax: number,
        options: {
            tolerance: number;
        }
    );

    eigenvalues(emin: number, emax: number): number[];

    eigenvaluesByIndex(imin: number, imax: number): number[];

    delete(): void;
}

declare interface matslise {
    Matslise: typeof Matslise;
    MatsliseHalfRange: typeof MatsliseHalfRange;
    Matslise2D: typeof Matslise2D;
}

export default class MatsliseModule {
    constructor(Module?: {
        locateFile(path: string, prefix: string): string;
    });

    then(resolve: (matslise: matslise) => void): void;
}
