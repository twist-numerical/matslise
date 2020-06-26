#ifndef MATSLISE_CHECKORTHONORMALITY_H
#define MATSLISE_CHECKORTHONORMALITY_H

#include "../catch.hpp"
#include "../../src/matslise.h"
#include "../../src/util/quadrature.h"

using namespace std;
using namespace Eigen;
using namespace matslise;
using namespace quadrature;

template<typename Scalar=double, typename doubleIterator>
void checkOrthonormality(const AbstractMatslise2D<Scalar> *p, const doubleIterator &begin, const doubleIterator &end) {
    const int n = 201;
    Array<Scalar, Dynamic, 1> x = lobatto::grid<Scalar>(
            Array<Scalar, Dynamic, 1>::LinSpaced(n, p->domain.getMin(0), p->domain.getMax(0)));
    Array<Scalar, Dynamic, 1> y = lobatto::grid<Scalar>(
            Array<Scalar, Dynamic, 1>::LinSpaced(n, p->domain.getMin(1), p->domain.getMax(1)));

    vector<Array<Scalar, Dynamic, Dynamic>> eigenfunctions;
    for (auto i = begin; i < end; ++i) {
        vector<Eigenfunction2D<>> fs = p->eigenfunction(*i);
        for (const Eigenfunction2D<> &f : fs) {
            ArrayXXd value = f(x, y);
            for (Index xi = 0; xi < n; xi += n / 20)
                for (Index yi = 0; yi < n; yi += n / 20)
                    REQUIRE(Approx(value(xi, yi)).margin(1e-10) == f(x[xi], y[yi]));
            eigenfunctions.push_back(value);
        }
    }

    auto eigenfunctionsBegin = eigenfunctions.begin();
    for (auto i = eigenfunctionsBegin; i != eigenfunctions.end(); ++i)
        for (auto j = eigenfunctionsBegin; j != eigenfunctions.end(); ++j) {
            INFO("Orthonormality of eigenfunction "
                         << std::distance(eigenfunctionsBegin, i) << " and " << std::distance(eigenfunctionsBegin, j));
            CHECK(Approx(lobatto::quadrature<Scalar>(x, y, *i * *j)).margin(1e-1) == (i == j ? 1 : 0));
        }

}


#endif //MATSLISE_CHECKORTHONORMALITY_H
