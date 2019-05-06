//
// Created by toon on 6/8/18.
//

#include <iostream>
#include "../matslise.h"
#include <memory>

using namespace std;
using namespace matslise;

#define EPS (1e-12)

HalfRange::HalfRange(function<double(double)> V, double xmax,
                     std::shared_ptr<matslise::SectorBuilder<Matslise>> sectorBuilder) {
    ms = new Matslise(V, 0, xmax, sectorBuilder);
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


inline bool isEven(const HalfRange *hr, double E, const matslise::Y<> &side, int even) {
    if (even == HalfRange::AUTO) {
        double error0 = get<0>(hr->ms->calculateError(E, Y<>({0, 1}, {0, 0}), side));
        double error1 = get<0>(hr->ms->calculateError(E, Y<>({0, 1}, {0, 0}), side));
        return abs(error0) < abs(error1);
    }
    return bool(even);
}

inline Y<> getY0(bool even) {
    return Y<>({even ? 1 : 0, even ? 0 : 1}, {0, 0});
}

Array<Y<>, Dynamic, 1>
HalfRange::computeEigenfunction(double E, const matslise::Y<> &side, const ArrayXd &x, int _even) const {
    long n = x.size();
    for (int i = 1; i < n; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");

    long negatives = 0;
    for (long i = 0; i < n; ++i)
        if (x[i] < 0)
            negatives = i + 1;

    ArrayXd xNeg(negatives);
    ArrayXd xPos(n - negatives);

    for (long i = 0; i < negatives; ++i)
        xNeg[i] = -x[negatives - 1 - i];

    for (long i = negatives; i < n; ++i)
        xPos[i - negatives] = x[i];


    bool even = isEven(this, E, side, _even);
    Y<> y = getY0(even);
    Array<Y<>, Dynamic, 1> yNeg, yPos, ys(n);
    yNeg = ms->computeEigenfunction(E, y, side, xNeg);
    yPos = ms->computeEigenfunction(E, y, side, xPos);

    for (long i = 0; i < negatives; ++i) {
        ys[negatives - 1 - i] = yNeg[i] * M_SQRT1_2;
        if (even)
            ys[negatives - 1 - i] *= -1;
    }
    for (long i = negatives; i < n; ++i)
        ys[i] = yPos[i - negatives] * M_SQRT1_2;

    return ys;
}


std::function<Y<>(double)> HalfRange::eigenfunctionCalculator(double E, const Y<> &side, int _even) const {
    bool even = isEven(this, E, side, _even);
    function<Y<>(double)> calculator(ms->eigenfunctionCalculator(E, getY0(even), side));
    return [calculator, even](double x) -> Y<> {
        Y<> c = calculator(x < 0 ? -x : x);
        c *= M_SQRT1_2;
        if (x < 0 && !even)
            c *= -1;
        return c;
    };
}

vector<pair<int, double>> *
mergeEigenvalues(vector<pair<int, double>> *even, vector<pair<int, double>> *odd) {
    vector<pair<int, double>> *values = new vector<pair<int, double>>();


    for (pair<int, double> &iE : *even)
        iE.first *= 2;

    for (pair<int, double> &iE : *odd)
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

double HalfRange::computeEigenvalueError(double E, const matslise::Y<> &side, int _even) const {
    return ms->computeEigenvalueError(E, getY0(isEven(this, E, side, _even)), side);
}

vector<pair<int, double>> *
HalfRange::computeEigenvaluesByIndex(int Imin, int Imax, const Y<> &side) const {
    return mergeEigenvalues(
            ms->computeEigenvaluesByIndex(Imin / 2 + Imin % 2, Imax / 2 + Imax % 2, Y<>({1, 0}, {0, 0}), side),
            ms->computeEigenvaluesByIndex(Imin / 2, Imax / 2, Y<>({0, 1}, {0, 0}), side));
}

vector<pair<int, double>> *
HalfRange::computeEigenvalues(double Emin, double Emax, const Y<> &side) const {
    return mergeEigenvalues(ms->computeEigenvalues(Emin, Emax, Y<>({1, 0}, {0, 0}), side),
                            ms->computeEigenvalues(Emin, Emax, Y<>({0, 1}, {0, 0}), side));
}

HalfRange::~HalfRange() {
    delete ms;
}
