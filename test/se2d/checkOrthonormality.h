#ifndef MATSLISE_CHECKORTHONORMALITY_H
#define MATSLISE_CHECKORTHONORMALITY_H

#include "../catch.hpp"
#include "../../src/matslise.h"
#include "../../src/util/lobatto.h"

using namespace std;
using namespace Eigen;
using namespace matslise;

template<typename Scalar, typename doubleIterator>
void checkOrthonormality(const SE2D<Scalar> &p, const doubleIterator &begin, const doubleIterator &end) {
    int n = 71;
    Array<Scalar, Dynamic, 1> x = lobatto::grid<Scalar>(
            Array<Scalar, Dynamic, 1>::LinSpaced(n, p.domain.getMin(0), p.domain.getMax(0)));
    Array<Scalar, Dynamic, 1> y = lobatto::grid<Scalar>(
            Array<Scalar, Dynamic, 1>::LinSpaced(n, p.domain.getMin(1), p.domain.getMax(1)));

    vector<Array<Scalar, Dynamic, Dynamic>> eigenfunctions;
    for (auto i = begin; i != end; ++i) {
        vector<Array<Scalar, Dynamic, Dynamic>> fs = p.computeEigenfunction(*i, x, y);
        for (const Array<Scalar, Dynamic, Dynamic> &f : fs)
            eigenfunctions.push_back(f);
    }

    for (auto i = eigenfunctions.begin(); i != eigenfunctions.end(); ++i)
        for (auto j = eigenfunctions.begin(); j != eigenfunctions.end(); ++j) {
            if (i != j) { // not needed when normalizing is done
                CHECKED_ELSE(Approx(lobatto::quadrature<Scalar>(x, y, *i * *j)).margin(1e-1) == (i == j ? 1 : 0)) {
                    auto l = begin;
                    for (auto k = eigenfunctions.begin(); k != eigenfunctions.end(); ++k) {
                        if (k == i)
                            break;
                        ++l;
                    }
                    FAIL(*l);
                }
            }
        }

}


#endif //MATSLISE_CHECKORTHONORMALITY_H
