#include <iostream>
#include "../matslise.h"
#include "../util/constants.h"
#include <memory>

using namespace std;
using namespace matslise;
using namespace Eigen;

#define EPS (1e-12)

template<typename Scaler>
HalfRange<Scaler>::HalfRange(
        function<Scaler(Scaler)> V, const Scalar &xmax, shared_ptr<SectorBuilder<Matslise<Scaler>>> sectorBuilder) {
    ms = new Matslise<Scaler>(V, 0, xmax, sectorBuilder);
}

template<typename T>
void removeDoubles(vector<T> &x) {
    sort(x.begin(), x.end());
    unsigned long s = x.size();
    int i = 0;
    for (unsigned long j = 1; j < s; ++j) {
        if (x[i] != x[j])
            x[++i] = x[j];
    }
    x.erase(x.begin() + i + 1, x.end());
}


template<typename Scalar>
inline bool isEven(const HalfRange<Scalar> *hr, Scalar E, const Y<Scalar> &side, int even) {
    if (even == HalfRange<Scalar>::AUTO) {
        Scalar error0 = get<0>(hr->ms->calculateError(E, Y<Scalar>({1, 0}, {0, 0}), side));
        Scalar error1 = get<0>(hr->ms->calculateError(E, Y<Scalar>({0, 1}, {0, 0}), side));
        return abs(error0) < abs(error1);
    }
    return bool(even);
}

template<typename Scalar>
inline Y<Scalar> getY0(bool even) {
    return Y<Scalar>({even ? 1 : 0, even ? 0 : 1}, {0, 0});
}

template<typename Scalar>
Array<Y<Scalar>, Dynamic, 1>
HalfRange<Scalar>::computeEigenfunction(const Scalar &E, const Y<Scalar> &side,
                                        const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x,
                                        int _even) const {
    long n = x.size();
    for (int i = 1; i < n; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");

    long negatives = 0;
    for (long i = 0; i < n; ++i)
        if (x[i] < 0)
            negatives = i + 1;

    Array<Scalar, Dynamic, 1> xNeg(negatives);
    Array<Scalar, Dynamic, 1> xPos(n - negatives);

    for (long i = 0; i < negatives; ++i)
        xNeg[i] = -x[negatives - 1 - i];

    for (long i = negatives; i < n; ++i)
        xPos[i - negatives] = x[i];


    bool even = isEven(this, E, side, _even);
    Y<Scalar> y = getY0<Scalar>(even);
    Array<Y<Scalar>, Dynamic, 1> yNeg, yPos, ys(n);
    yNeg = ms->computeEigenfunction(E, y, side, xNeg);
    yPos = ms->computeEigenfunction(E, y, side, xPos);

    static const Scalar SQRT1_2 = sqrt(Scalar(.5));
    for (long i = 0; i < negatives; ++i) {
        Scalar f = (even ? 1 : -1) * SQRT1_2;
        ys[negatives - 1 - i].y = f * DiagonalMatrix<Scalar, 2>(1, -1) * yNeg[i].y;
        ys[negatives - 1 - i].dy = f * DiagonalMatrix<Scalar, 2>(1, -1) * yNeg[i].dy;
    }
    for (long i = negatives; i < n; ++i)
        ys[i] = yPos[i - negatives] * SQRT1_2;

    return ys;
}


template<typename Scalar>
std::function<Y<Scalar>(Scalar)>
HalfRange<Scalar>::eigenfunctionCalculator(const Scalar &E, const Y<Scalar> &side, int _even) const {
    bool even = isEven(this, E, side, _even);
    function<Y<Scalar>(Scalar)> calculator(ms->eigenfunctionCalculator(E, getY0<Scalar>(even), side));
    return [calculator, even](Scalar x) -> Y<Scalar> {
        Y<Scalar> c = calculator(x < 0 ? -x : x);
        c *= M_SQRT1_2;
        if (x < 0 && !even)
            c *= -1;
        return c;
    };
}

template<typename Scalar>
vector<pair<int, Scalar>> *
mergeEigenvalues(vector<pair<int, Scalar>> *even, vector<pair<int, Scalar>> *odd) {
    vector<pair<int, Scalar>> *values = new vector<pair<int, Scalar>>();


    for (pair<int, Scalar> &iE : *even)
        iE.first *= 2;

    for (pair<int, Scalar> &iE : *odd)
        iE.first = 2 * iE.first + 1;

    auto a = even->begin();
    auto b = odd->begin();
    while (a != even->end() || b != odd->end()) {
        if (a == even->end())
            values->push_back(*b++);
        else if (b == odd->end())
            values->push_back(*a++);
        else {
            if (*a < *b)
                values->push_back(*a++);
            else
                values->push_back(*b++);
        }
    }

    delete even;
    delete odd;
    return values;

}

template<typename Scalar>
Scalar HalfRange<Scalar>::computeEigenvalueError(const Scalar &E, const Y<Scalar> &side, int _even) const {
    return ms->computeEigenvalueError(E, getY0<Scalar>(isEven(this, E, side, _even)), side);
}

template<typename Scalar>
vector<pair<int, Scalar>> *
HalfRange<Scalar>::computeEigenvaluesByIndex(int Imin, int Imax, const Y<Scalar> &side, const Y<Scalar> &right) const {
    checkSymmetry(side, right);
    return mergeEigenvalues(
            ms->computeEigenvaluesByIndex(Imin / 2 + Imin % 2, Imax / 2 + Imax % 2, Y<Scalar>({1, 0}, {0, 0}), side),
            ms->computeEigenvaluesByIndex(Imin / 2, Imax / 2, Y<Scalar>({0, 1}, {0, 0}), side));
}

template<typename Scalar>
vector<pair<int, Scalar>> *
HalfRange<Scalar>::computeEigenvalues(
        const Scalar &Emin, const Scalar &Emax, const Y<Scalar> &side, const Y<Scalar> &right) const {
    checkSymmetry(side, right);
    return mergeEigenvalues(ms->computeEigenvalues(Emin, Emax, Y<Scalar>({1, 0}, {0, 0}), side),
                            ms->computeEigenvalues(Emin, Emax, Y<Scalar>({0, 1}, {0, 0}), side));
}

template<typename Scalar>
HalfRange<Scalar>::~HalfRange() {
    delete ms;
}

#include "../util/instantiate.h"