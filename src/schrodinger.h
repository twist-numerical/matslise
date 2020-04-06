#ifndef MATSLISE_SCHRODINGER_H
#define MATSLISE_SCHRODINGER_H

#include "util/sectorbuilder.h"

namespace matslise {
    template<typename Scalar=double>
    class Matslise2D;

    template<typename Scalar=double>
    struct Options1 {
        SectorBuilder<matslise::Matslise<Scalar>, Scalar> _builder = matslise::sector_builder::uniform<Matslise<Scalar>>(26);
        bool _symmetric = false;

        Options1<Scalar> &symmetric(bool symmetric) {
            _symmetric = symmetric;
            return *this;
        }

        Options1<Scalar> &sectorCount(int count) {
            _builder = matslise::sector_builder::uniform<Matslise<Scalar>>(count);
            return *this;
        }

        Options1<Scalar> &tolerance(Scalar tol) {
            _builder = matslise::sector_builder::automatic<Matslise<Scalar>>(tol);
            return *this;
        }
    };

    template<typename Scalar=double>
    struct Options2 {
        Options1<Scalar> nestedOptions;
        matslise::SectorBuilder<matslise::Matslise2D<Scalar>, Scalar> _builder = matslise::sector_builder::uniform<Matslise2D<Scalar>>(17);
        int _stepsPerSector = 1;
        int _N = 12;
        int _gridPoints = 52;

        Options2<Scalar> &nested(Options1<Scalar> &nested) {
            nestedOptions = nested;
            return *this;
        }

        Options2<Scalar> &sectorCount(int count) {
            _builder = matslise::sector_builder::uniform<Matslise2D<Scalar>>(count);
            return *this;
        }

        Options2<Scalar> &tolerance(Scalar tol) {
            _builder = matslise::sector_builder::automatic<Matslise2D<Scalar>>(tol);
            return *this;
        }

        Options2<Scalar> &stepsPerSector(int count) {
            _stepsPerSector = count;
            return *this;
        }

        Options2<Scalar> &N(int value) {
            _N = value;
            return *this;
        }

        Options2<Scalar> &gridPoints(int value) {
            _gridPoints = value;
            return *this;
        }
    };


}

#endif //MATSLISE_SCHRODINGER_H
