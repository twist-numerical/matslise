#include <iostream>
#include "../matslise.h"
#include "../util/constants.h"
#include <memory>

using namespace std;
using namespace matslise;
using namespace Eigen;

#define EPS (1e-12)

template<typename Scalar>
MatsliseHalf<Scalar>::MatsliseHalf(function<Scalar(Scalar)> V, const Scalar &xmax, const Scalar &tolerance,
                                   SectorBuilder<Matslise<Scalar>> sectorBuilder) {
    ms = new Matslise<Scalar>(V, 0, xmax, tolerance, sectorBuilder);
}

template<typename Scalar>
inline bool isEven(const MatsliseHalf<Scalar> *hr, Scalar E, const Y<Scalar> &side, int index) {
    if (index == -1) {
        Scalar error0 = get<0>(hr->ms->matchingError(E, Y<Scalar>({1, 0}, {0, 0}), side));
        Scalar error1 = get<0>(hr->ms->matchingError(E, Y<Scalar>({0, 1}, {0, 0}), side));
        return abs(error0) < abs(error1);
    }
    return bool(index % 2 == 0);
}

template<typename Scalar>
inline Y<Scalar> getY0(bool even) {
    return Y<Scalar>({even ? 1 : 0, even ? 0 : 1}, {0, 0});
}

template<typename Scalar>
Array<Y<Scalar>, Dynamic, 1>
MatsliseHalf<Scalar>::eigenfunction(
        const Scalar &E, const Y<Scalar> &side, const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x, int index) const {
    Eigen::Index n = x.size();
    for (Eigen::Index i = 1; i < n; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");

    Eigen::Index negatives = 0;
    for (Eigen::Index i = 0; i < n; ++i)
        if (x[i] < 0)
            negatives = i + 1;

    Array<Scalar, Dynamic, 1> xNeg(negatives);
    Array<Scalar, Dynamic, 1> xPos(n - negatives);

    for (Eigen::Index i = 0; i < negatives; ++i)
        xNeg[i] = -x[negatives - 1 - i];

    for (Eigen::Index i = negatives; i < n; ++i)
        xPos[i - negatives] = x[i];


    bool even = isEven(this, E, side, index);
    Y<Scalar> y = getY0<Scalar>(even);
    Array<Y<Scalar>, Dynamic, 1> yNeg, yPos, ys(n);
    yNeg = ms->eigenfunction(E, y, side, xNeg);
    yPos = ms->eigenfunction(E, y, side, xPos);

    static const Scalar SQRT1_2 = sqrt(Scalar(.5));
    for (Eigen::Index i = 0; i < negatives; ++i) {
        Scalar f = (even ? 1 : -1) * SQRT1_2;
        ys[negatives - 1 - i].y = f * DiagonalMatrix<Scalar, 2>(1, -1) * yNeg[i].y;
        ys[negatives - 1 - i].dy = f * DiagonalMatrix<Scalar, 2>(1, -1) * yNeg[i].dy;
    }
    for (Eigen::Index i = negatives; i < n; ++i)
        ys[i] = yPos[i - negatives] * SQRT1_2;

    return ys;
}


template<typename Scalar>
std::function<Y<Scalar>(Scalar)>
MatsliseHalf<Scalar>::eigenfunctionCalculator(const Scalar &E, const Y<Scalar> &side, int index) const {
    bool even = isEven(this, E, side, index);
    function<Y<Scalar>(Scalar)> calculator(ms->eigenfunctionCalculator(E, getY0<Scalar>(even), side));
    return [calculator, even](Scalar x) -> Y<Scalar> {
        Y<Scalar> c = calculator(x < 0 ? -x : x);
        c *= sqrt(Scalar(.5));
        if (x < 0 && !even)
            c *= -1;
        return c;
    };
}

template<typename Scalar>
vector<pair<int, Scalar>>
mergeEigenvalues(const vector<pair<int, Scalar>> &even, const vector<pair<int, Scalar>> &odd) {
    vector<pair<int, Scalar>> values;

    auto a = even.begin();
    auto b = odd.begin();
    while (a != even.end() || b != odd.end()) {
        if (a == even.end() || (b != odd.end() && b->first < a->first)) {
            values.push_back({2 * b->first + 1, b->second});
            ++b;
        } else {
            values.push_back({2 * a->first, a->second});
            ++a;
        }
    }

    return values;
}

template<typename Scalar>
Scalar MatsliseHalf<Scalar>::eigenvalueError(const Scalar &E, const Y<Scalar> &side, int index) const {
    return ms->eigenvalueError(E, getY0<Scalar>(isEven(this, E, side, index)), side);
}

template<typename Scalar>
vector<pair<int, Scalar>>
MatsliseHalf<Scalar>::eigenvaluesByIndex(int Imin, int Imax, const Y<Scalar> &side) const {
    return mergeEigenvalues(
            ms->eigenvaluesByIndex(Imin / 2 + Imin % 2, Imax / 2 + Imax % 2, Y<Scalar>({1, 0}, {0, 0}), side),
            ms->eigenvaluesByIndex(Imin / 2, Imax / 2, Y<Scalar>({0, 1}, {0, 0}), side));
}

template<typename Scalar>
vector<pair<int, Scalar>>
MatsliseHalf<Scalar>::eigenvalues(
        const Scalar &Emin, const Scalar &Emax, const Y<Scalar> &side) const {
    return mergeEigenvalues(ms->eigenvalues(Emin, Emax, Y<Scalar>({1, 0}, {0, 0}), side),
                            ms->eigenvalues(Emin, Emax, Y<Scalar>({0, 1}, {0, 0}), side));
}

template<typename Scalar>
MatsliseHalf<Scalar>::~MatsliseHalf() {
    delete ms;
}

#include "../util/instantiate.h"