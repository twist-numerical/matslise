#ifndef MATSLISE_CHECKORTHONORMALITY_H
#define MATSLISE_CHECKORTHONORMALITY_H

#include "../catch.hpp"
#include "../../matslise/matslise.h"
#include "../../matslise/util/quadrature.h"

using namespace std;
using namespace Eigen;
using namespace matslise;
using namespace quadrature;

template<typename Scalar=double>
void checkOrthonormality(const AbstractMatslise2D<Scalar> &p, const vector<Eigenfunction2D<Scalar>> &eigenfunctions) {
    const int n = 101;
    Array<Scalar, Dynamic, 1> x = lobatto::grid<Scalar>(
            Array<Scalar, Dynamic, 1>::LinSpaced(n, p.domain.min(0), p.domain.max(0)));
    Array<Scalar, Dynamic, 1> y = lobatto::grid<Scalar>(
            Array<Scalar, Dynamic, 1>::LinSpaced(n, p.domain.min(1), p.domain.max(1)));

    vector<Array<Scalar, Dynamic, Dynamic>> eigenfunctionsOnGrid;
    for (unsigned long i = 0; i < eigenfunctions.size(); ++i) {
        const Eigenfunction2D<Scalar> &f = eigenfunctions[i];
        ArrayXXd value = f(x, y);
        for (Index xi = 0; xi < n; xi += n / 20)
            for (Index yi = 0; yi < n; yi += n / 20)
                REQUIRE(Approx(value(xi, yi)).margin(1e-10) == f(x[xi], y[yi]));
        eigenfunctionsOnGrid.push_back(value);
    }

    auto eigenfunctionsBegin = eigenfunctionsOnGrid.begin();
    for (auto i = eigenfunctionsBegin; i != eigenfunctionsOnGrid.end(); ++i)
        for (auto j = eigenfunctionsBegin; j != eigenfunctionsOnGrid.end(); ++j) {
            INFO("Orthonormality of eigenfunction "
                         << std::distance(eigenfunctionsBegin, i) << " and " << std::distance(eigenfunctionsBegin, j));
            CHECK(Approx(lobatto::quadrature<Scalar>(x, y, *i * *j)).margin(1e-2) == (i == j ? 1 : 0));
        }

}

template<typename Scalar=double>
void checkProblem(const AbstractMatslise2D<Scalar> &p, const vector<tuple<Index, Scalar, Index>> &exactEigenvalues,
                  const Scalar &tolerance = 1e-7) {
    Index count = get<0>(exactEigenvalues.back()) + get<2>(exactEigenvalues.back());
    Index first = get<0>(exactEigenvalues.front());

    vector<tuple<Index, Scalar, Index>> eigenvalues = p.eigenvaluesByIndex(first, first + count);

    auto exact = exactEigenvalues.begin();
    auto found = eigenvalues.begin();
    vector<Eigenfunction2D<Scalar>> eigenfunctions;
    for (; exact != exactEigenvalues.end() && found != eigenvalues.end(); ++exact, ++found) {
        INFO("Checking eigenvalue: (" << get<0>(*exact) << ", " << get<1>(*exact) << ", " << get<2>(*exact)
                                      << ") against (" << get<0>(*found) << ", " << get<1>(*found) << ", "
                                      << get<2>(*found) << ")")
        CHECK(get<0>(*exact) == get<0>(*found));
        CHECK(Approx(get<1>(*exact)).margin(tolerance) == get<1>(*found));
        CHECK(get<2>(*exact) == get<2>(*found));
        auto currentEigenfunctions = p.eigenfunction(get<1>(*found));
        CHECK(get<2>(*exact) == (Index) currentEigenfunctions.size());
        for (auto &eigenfunction : currentEigenfunctions)
            eigenfunctions.push_back(eigenfunction);
    }

    checkOrthonormality(p, eigenfunctions);
}


#endif //MATSLISE_CHECKORTHONORMALITY_H
