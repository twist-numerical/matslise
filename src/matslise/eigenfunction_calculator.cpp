//
// Created by toon on 9/11/18.
//

#include <functional>
#include "../matslise.h"

using namespace matslise;
using namespace matslise::matslise_util;

EigenfunctionCalculator::EigenfunctionCalculator(const Matslise *ms, double E, const Y<> &left, const Y<> &right)
        : ms(ms), E(E) {
    ys.reserve(ms->sectorCount + 1);
    int m = ms->sectorCount / 2 + 1;
    ys[0] = left;
    for (int i = 1; i <= m; ++i)
        ys[i] = ms->sectors[i - 1]->propagate(E, ys[i - 1], true);
    ys[ms->sectorCount] = right;
    for (int i = ms->sectorCount - 1; i > m; --i)
        ys[i] = ms->sectors[i]->propagate(E, ys[i + 1], false);
    Y<> yr = ms->sectors[m]->propagate(E, ys[m + 1], false);
    double s = ys[m].y[0] / yr.y[0];
    for (int i = m + 1; i < ms->sectorCount; ++i)
        ys[i] *= s;
}

Y<> EigenfunctionCalculator::operator()(double x) const {
    int a = 0;
    int b = ms->sectorCount;
    while (a + 1 < b) {
        int c = (a + b) / 2;
        if (x < ms->sectors[c]->xmin)
            b = c;
        else
            a = c;
    }
    return ms->sectors[a]->propagate(E, ys[a], x - ms->sectors[a]->xmin);
}

EigenfunctionCalculator::~EigenfunctionCalculator() {
}
