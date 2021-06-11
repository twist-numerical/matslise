#ifndef MATSLISE_THETA_H
#define MATSLISE_THETA_H

#include "../matslise.h"
#include "constants.h"

template<typename Scalar>
inline Scalar atan_safe(const Scalar &y, const Scalar &x) {
    Scalar r = atan2(y, x);
    if (r > constants<Scalar>::PI / 2)
        r -= constants<Scalar>::PI;
    else if (r <= -constants<Scalar>::PI / 2)
        r += constants<Scalar>::PI;
    return r;
}

namespace matslise {
    template<typename Scalar>
    inline Scalar theta(const Y <Scalar> &y) {
        return atan_safe<Scalar>(y.y(0, 0), y.y(1, 0));
    }
}

#endif