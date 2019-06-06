//
// Created by toon on 6/6/19.
//

#ifndef MATSLISE_MATH_H
#define MATSLISE_MATH_H

#include <cmath>

#define FMATH_FUNCTION(name, function) \
static Scalar name(const Scalar &a) {\
    return function(a);\
}

#define FMATH_SPECIAL(type, name, function) \
namespace std {\
    type name(const float128 &a) {\
        return function(a);\
    }\
}\
template<> type fmath<float128>::name(const float128 &a) {\
    return function(a);\
}

template<typename Scalar>
class fmath {
public:
    static const Scalar PI;

    FMATH_FUNCTION(abs, std::fabs)

    FMATH_FUNCTION(sqrt, std::sqrt)

    FMATH_FUNCTION(sin, std::sin)

    FMATH_FUNCTION(cos, std::cos)

    FMATH_FUNCTION(sinh, std::sinh)

    FMATH_FUNCTION(cosh, std::cosh)

    FMATH_FUNCTION(atan, std::atan)

    static Scalar pow(const Scalar &a, const Scalar &b) {
        return std::pow(a, b);
    }

    static bool isnan(const Scalar &a) {
        return std::isnan((double) a);
    }

    static int ceil(const Scalar &a) {
        return (int) std::ceil((double) a);
    }

    static int floor(const Scalar &a) {
        return (int) std::floor((double) a);
    }

    static int round(const Scalar &a) {
        return (int) std::round((double) a);
    }
};

#endif