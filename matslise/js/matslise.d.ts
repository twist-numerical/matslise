type Eigenfunction = {
    (x: number): [number, number];

    (x: number[]): [number, number][];

    delete(): void;
};

declare class AbstractMatslise {
    eigenvaluesByIndex(
        imin: number,
        imax: number,
        left: [number, number],
        right: [number, number]
    ): { index: number; eigenvalue: number }[];

    eigenpairsByIndex(
        imin: number,
        imax: number,
        left: [number, number],
        right: [number, number]
    ): { index: number; eigenvalue: number; eigenfunction: Eigenfunction }[];

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

declare interface matslise {
    Matslise: typeof Matslise;
    MatsliseHalf: typeof MatsliseHalf;
}

export default class MatsliseModule {
    constructor(Module?: {
        locateFile(path: string, prefix: string): string;
    });

    then(resolve: (matslise: matslise) => void): void;
}
