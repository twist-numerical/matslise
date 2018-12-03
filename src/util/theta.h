#ifndef MATSLISE_THETA_H
#define MATSLISE_THETA_H

#include "../matslise.h"

inline double atan_safe(double y, double x) {
    if (x == 0)
        return y == 0 ? 0 : M_PI_2;
    return atan(y / x);
}

namespace matslise {
    inline double theta(const Y<double> &y) {
        return atan_safe(y.y[0], y.y[1]);
    }
}

#endif