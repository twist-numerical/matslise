//
// Created by toon on 6/8/18.
//

#include <iostream>
#include "../matslise.h"

using namespace std;
using namespace matslise;

#define EPS (1e-12)

HalfRange::HalfRange(function<double(double)> V, double xmax, int sectorCount) {
    ms = new Matslise(V, 0, xmax, sectorCount);
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

Array<Y<double>, Dynamic, 1>
HalfRange::computeEigenfunction(double E, const matslise::Y<double> &side, const ArrayXd &x) const {
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


    Y<double> y0({0, 1});
    Y<double> y1({1, 0});

    double error0 = get<0>(ms->calculateError(E, y0, side));
    double error1 = get<0>(ms->calculateError(E, y1, side));

    Array<Y<double>, Dynamic, 1> yNeg, yPos, ys(n);
    bool is0 = abs(error0) < abs(error1);

    yNeg = ms->computeEigenfunction(E, is0 ? y0 : y1, side, xNeg);
    yPos = ms->computeEigenfunction(E, is0 ? y0 : y1, side, xPos);

    for (long i = 0; i < negatives; ++i)
        ys[negatives - 1 - i] = is0 ? -yNeg[i] : yNeg[i];
    for (long i = negatives; i < n; ++i)
        ys[i] = yPos[i - negatives];

    return ys;
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

};

vector<pair<int, double>> *
HalfRange::computeEigenvaluesByIndex(int Imin, int Imax, const Y<double> &side) const {
    return mergeEigenvalues(
            ms->computeEigenvaluesByIndex(Imin / 2 + Imin % 2, Imax / 2 + Imax % 2, Y<double>({1, 0}), side),
            ms->computeEigenvaluesByIndex(Imin / 2, Imax / 2, Y<double>({0, 1}), side));
};

vector<pair<int, double>> *
HalfRange::computeEigenvalues(double Emin, double Emax, const Y<double> &side) const {
    return mergeEigenvalues(ms->computeEigenvalues(Emin, Emax, Y<double>({1, 0}), side),
                            ms->computeEigenvalues(Emin, Emax, Y<double>({0, 1}), side));
}

HalfRange::~HalfRange() {
    delete ms;
}
