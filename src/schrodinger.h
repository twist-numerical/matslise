//
// Created by toon on 12/7/18.
//

#ifndef MATSLISE_SCHRODINGER_H
#define MATSLISE_SCHRODINGER_H

#include "matslise.h"
#include "se2d.h"

using namespace Eigen;
using namespace std;

namespace matslise {

    template<int n>
    struct Options {
        Options<n - 1> nestedOptions;
        int _sectorCount = 17;
        int _stepsPerSector = 4;
        int _N = 12;
        int _gridPoints = 52;

        Options<n> &nested(Options<n - 1> &nested) {
            nestedOptions = nested;
            return *this;
        }

        Options<n> &sectorCount(int count) {
            _sectorCount = count;
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
        int _sectorCount = 26;

        Options<1> &sectorCount(int count) {
            _sectorCount = count;
            return *this;
        }
    };

    template<int n>
    class Rectangle {
    public:
        Rectangle<n - 1> sub;
        double min, max;

        double getMin(int axis) const {
            if(axis+1 == n)
                return min;
            return sub.getMin(axis);
        }

        double getMax(int axis) const {
            if(axis+1 == n)
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
            (void)(axis); // UNUSED
            return min;
        }

        double getMax(int axis) const {
            assert(axis == 0);
            (void)(axis); // UNUSED
            return max;
        }
    };
}

#endif //MATSLISE_SCHRODINGER_H
