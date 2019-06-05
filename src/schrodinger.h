//
// Created by toon on 12/7/18.
//

#ifndef MATSLISE_SCHRODINGER_H
#define MATSLISE_SCHRODINGER_H

#include "matslise.h"
#include "se2d.h"

using namespace Eigen;

namespace matslise {
    template<int n>
    class SEnD;

    template<int n>
    struct Options {
        Options<n - 1> nestedOptions;
        std::shared_ptr<matslise::SectorBuilder<matslise::SEnD<n>>> _builder
                = std::shared_ptr<matslise::SectorBuilder<matslise::SEnD<n>>>(
                        new matslise::sectorbuilder::Uniform<SEnD<n>>(17));
        int _stepsPerSector = 1;
        int _N = 12;
        int _gridPoints = 52;

        Options<n> &nested(Options<n - 1> &nested) {
            nestedOptions = nested;
            return *this;
        }

        Options<n> &sectorCount(int count) {
            _builder.reset(new matslise::sectorbuilder::Uniform<SEnD<n>>(count));
            return *this;
        }

        Options<n> &tolerance(double tol) {
            _builder.reset(new matslise::sectorbuilder::Auto<SEnD<n>>(tol));
            return *this;
        }

        Options<n> &stepsPerSector(int count) {
            _stepsPerSector = count;
            return *this;
        }

        Options<n> &N(int value) {
            _N = value;
            return *this;
        }

        Options<n> &gridPoints(int value) {
            _gridPoints = value;
            return *this;
        }
    };

    template<>
    struct Options<1> {
        std::shared_ptr<matslise::SectorBuilder<matslise::Matslise>> _builder = Matslise::UNIFORM(26);
        bool _symmetric = false;

        Options<1> &symmetric(bool symmetric) {
            _symmetric = symmetric;
            return *this;
        }

        Options<1> &sectorCount(int count) {
            _builder = Matslise::UNIFORM(count);
            return *this;
        }

        Options<1> &tolerance(double tol) {
            _builder = Matslise::AUTO(tol);
            return *this;
        }
    };

    template<int n>
    class Rectangle {
    public:
        Rectangle<n - 1> sub;
        double min, max;

        double getMin(int axis) const {
            if (axis + 1 == n)
                return min;
            return sub.getMin(axis);
        }

        double getMax(int axis) const {
            if (axis + 1 == n)
                return max;
            return sub.getMax(axis);
        }
    };

    template<>
    class Rectangle<1> {
    public:
        double min, max;

        double getMin(int axis) const {
            assert(axis == 0);
            (void) (axis); // UNUSED
            return min;
        }

        double getMax(int axis) const {
            assert(axis == 0);
            (void) (axis); // UNUSED
            return max;
        }
    };
}

#endif //MATSLISE_SCHRODINGER_H
