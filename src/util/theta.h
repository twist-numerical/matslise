#ifndef MATSLISE_THETA_H
#define MATSLISE_THETA_H

#include "../matslise.h"

inline double atan_pos(double y, double x) {
    if (x == 0)
        return y == 0 ? 0 : M_PI_2;
    double r = atan(y / x);
    if (r < 0)
        r += M_PI;
    return r;
}

namespace matslise {
    inline double theta(const Y<double> &y) {
        return atan_pos(y.y[0], y.y[1]);
    }
}

#endif