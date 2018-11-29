//
// Created by toon on 9/11/18.
//

#include "../matscs.h"
#include <array>

using namespace matslise;
using namespace matslise::matscs_util;

EigenfunctionCalculator::EigenfunctionCalculator(Matscs *ms, double E, const Y<VectorXd> &left, const Y<VectorXd> &right)
        : ms(ms), E(E) {
    ys.reserve(ms->sectorCount + 1);
    int m = ms->sectorCount / 2 + 1;
    ys[0] = left;
    for (int i = 1; i <= m; ++i)
        ys[i] = ms->sectors[i - 1]->propagate(E, ys[i - 1], true);
    ys[ms->sectorCount] = right;
    for (int i = ms->sectorCount - 1; i > m; --i)
        ys[i] = ms->sectors[i]->propagate(E, ys[i + 1], false);
    Y<VectorXd> yr = ms->sectors[m]->propagate(E, ys[m + 1], false);
    double s = ys[m].y[0].norm() / yr.y[0].norm();
    for (int i = m + 1; i < ms->sectorCount; ++i)
        ys[i] *= s;
}

Y<VectorXd> EigenfunctionCalculator::eval(double x) const {
    int a = 0;
    int b = ms->sectorCount;
    while (a + 1 < b) {
        int c = (a + b) / 2;
        if (x < ms->sectors[c]->xmin)
            b = c;
        else
            a = c;
    }
    return ms->sectors[a]->propagate<double, -1, 1>(E, ys[a], x - ms->sectors[a]->xmin);
}
