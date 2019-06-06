#ifndef MATSLISE_THETA_H
#define MATSLISE_THETA_H

#include "../matslise.h"
#include "fmath.h"

template<typename Scalar>
inline Scalar atan_safe(Scalar y, Scalar x) {
    if (x == 0)
        return y == 0 ? 0 : fmath<Scalar>::PI / 2;
    return fmath<Scalar>::atan(y / x);
}

namespace matslise {
    template<typename Scalar>
    inline Scalar theta(const Y <Scalar> &y) {
        return atan_safe<Scalar>(y.y(0, 0), y.y(1, 0));
    }
}

#endif