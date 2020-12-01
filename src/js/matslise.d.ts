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
        tolerance: number
    );
}

declare class MatsliseHalf extends AbstractMatslise {
    constructor(
        potential: (x: number) => number,
        xmax: number,
        tolerance: number
    );
}

declare class AbstractMatslise2D {

    eigenvalue(E: number): [number, number];

    eigenvalues(emin: number, emax: number): { index: number, value: number, multiplicity: number }[];

    eigenvaluesByIndex(imin: number, imax: number): { index: number, value: number, multiplicity: number }[];

    eigenvalueError(E: number): number;

    computeEigenfunction(E: number, x: number[], y: number[]): number[][][];

    delete(): void;
}

declare class Matslise2D extends AbstractMatslise2D {
    constructor(
        potential: (x: number, y: number) => number,
        xmin: number,
        xmax: number,
        ymin: number,
        ymax: number,
        options: {
            tolerance?: number;
            xSectorCount?: number;
            xTolerance?: number;
            ySectorCount?: number;
            yTolerance?: number;
            xSymmetric?: boolean;
        }
    );
}

declare class Matslise2DHalf extends AbstractMatslise2D {
    constructor(
        potential: (x: number, y: number) => number,
        xmin: number,
        xmax: number,
        ymax: number,
        options: {
            tolerance?: number;
            xSectorCount?: number;
            xTolerance?: number;
            ySectorCount?: number;
            yTolerance?: number;
            xSymmetric?: boolean;
        }
    );
}

declare interface matslise {
    Matslise: typeof Matslise;
    MatsliseHalf: typeof MatsliseHalf;
    Matslise2D: typeof Matslise2D;
    Matslise2DHalf: typeof Matslise2DHalf;
}

export default class MatsliseModule {
    constructor(Module?: {
        locateFile(path: string, prefix: string): string;
    });

    then(resolve: (matslise: matslise) => void): void;
}
