//
// Created by toon on 6/8/18.
//

#include <iostream>
#include <algorithm>
#include "../matslise.h"

using namespace std;
using namespace matslise;

HalfRange::HalfRange(function<double(double)> V, double xmax, int sectorCount) {
    ms = new Matslise([V](double x) -> double {
        if (x < 0)
            return V(-x);
        else
            return V(x);
    }, 0, xmax, sectorCount);
}

template <typename T>
void removeDoubles(vector<T> &x) {
    sort(x.begin(), x.end());
    unsigned long s = x.size();
    int i = 0;
    for(unsigned long j = 1; j < s; ++j) {
        if(x[i] != x[j])
            x[++i] = x[j];
    }
    x.erase(x.begin()+i+1, x.end());
}

vector<Y> *HalfRange::computeEigenfunction(double E, const Y &side, vector<double> &x) const {
    for (auto i = x.begin(); i != x.end(); ++i)
        *i = abs(*i);
    removeDoubles(x);

    Y y0({0, 1});
    Y y1({1, 0});

    double error0 = get<0>(ms->calculateError(E, y0, side));
    double error1 = get<0>(ms->calculateError(E, y1, side));
    vector<Y> * y;
    bool is0 = abs(error0) < abs(error1);
    if(is0) {
        y = ms->computeEigenfunction(E, y0, side, x);
    } else {
        y = ms->computeEigenfunction(E, y1, side, x);
    }

    unsigned long xs = x.size();
    x.insert(x.begin(), x.rbegin(), x.rend());
    for(unsigned long i = 0; i < xs; ++i)
        x[i] = -x[i];
    y->insert(y->begin(), y->rbegin(), y->rend());

    if(is0)
        for(unsigned long i = 0; i < xs; ++i) {
            (*y)[i].y *= -1;
            (*y)[i].dy *= -1;
        }

    return y;
}

vector<double> *HalfRange::computeEigenvalues(double Emin, double Emax, const Y &side) const {
    vector<double> *even = ms->computeEigenvalues(Emin, Emax, Y({1, 0}), side);
    vector<double> *odd = ms->computeEigenvalues(Emin, Emax, Y({0, 1}), side);

    vector<double> *values = new vector<double>();

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

HalfRange::~HalfRange() {
    delete ms;
}
