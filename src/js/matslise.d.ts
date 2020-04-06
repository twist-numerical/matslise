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

    firstEigenvalue(): number;

    eigenvalue(E: number): number;

    eigenvalues(emin: number, emax: number): number[];

    eigenvaluesByIndex(imin: number, imax: number): number[];

    eigenvalueError(E: number): number;

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
            sectorCount?: number;
            nested?: {
                symmetric?: boolean;
                tolerance?: number;
                sectorCount?: number;
            }
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
            sectorCount?: number;
            nested?: {
                symmetric?: boolean;
                tolerance?: number;
                sectorCount?: number;
            }
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
