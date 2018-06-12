//
// Created by toon on 6/8/18.
//

#include <iostream>
#include <algorithm>
#include "../matslise.h"

using namespace std;
using namespace matslise;

#define EPS (1e-12)

HalfRange::HalfRange(function<double(double)> V, double xmax, int sectorCount) {
    ms = new Matslise([V](double x) -> double {
        if (x < 0)
            return V(-x);
        else
            return V(x);
    }, 0, xmax, sectorCount);
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

vector<Y> *HalfRange::computeEigenfunction(double E, const Y &side, vector<double> &x) const {
    sort(x.begin(), x.end());

    vector<double> absx = x;
    for (auto i = absx.begin(); i != absx.end(); ++i)
        *i = abs(*i);
    removeDoubles(absx);

    Y y0({0, 1});
    Y y1({1, 0});

    double error0 = get<0>(ms->calculateError(E, y0, side));
    double error1 = get<0>(ms->calculateError(E, y1, side));
    vector<Y> *absy;
    bool is0 = abs(error0) < abs(error1);
    if (is0)
        absy = ms->computeEigenfunction(E, y0, side, absx);
    else
        absy = ms->computeEigenfunction(E, y1, side, absx);


    vector<Y> *ys = new vector<Y>();
    auto i = x.begin();
    {
        auto absi = absx.rbegin();
        auto absj = absy->rbegin();
        while (i != x.end() && *i < 0) {
            while (*i > -*absi) {
                absj++;
                absi++;
            }
            if (*i == -*absi) {
                ys->push_back(is0 ? -*absj : *absj);
                ++i;
            } else
                i = x.erase(i);
        }
    }
    {
        auto absi = absx.begin();
        auto absj = absy->begin();
        while (i != x.end()) {
            while (*i > *absi) {
                absj++;
                absi++;
            }
            if (*i == *absi) {
                ys->push_back(*absj);
                ++i;
            } else
                i = x.erase(i);
        }
    }

    delete absy;

    return ys;
}

vector<tuple<unsigned int, double>> *mergeEigenvalues(vector<tuple<unsigned int, double>> *even, vector<tuple<unsigned int, double>> *odd) {
    vector<tuple<unsigned int, double>> *values = new vector<tuple<unsigned int, double>>();


    for (tuple<unsigned int, double> &iE : *even)
        get<0>(iE) *= 2;

    for (tuple<unsigned int, double> &iE : *odd)
        get<0>(iE) = 2 * get<0>(iE) + 1;

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

vector<tuple<unsigned int, double>> *HalfRange::computeEigenvaluesByIndex(unsigned int Imin, unsigned int Imax, const Y &side) const {
    return mergeEigenvalues(
            ms->computeEigenvaluesByIndex(Imin/2 + Imin%2, Imax/2 + Imax%2, Y({1, 0}), side),
            ms->computeEigenvaluesByIndex(Imin/2, Imax/2, Y({0, 1}), side));
};

vector<tuple<unsigned int, double>> *HalfRange::computeEigenvalues(double Emin, double Emax, const Y &side) const {
    return mergeEigenvalues(ms->computeEigenvalues(Emin, Emax, Y({1, 0}), side), ms->computeEigenvalues(Emin, Emax, Y({0, 1}), side));
}

HalfRange::~HalfRange() {
    delete ms;
}
