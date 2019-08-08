#ifndef MATSLISE_CONSTANTS_H
#define MATSLISE_CONSTANTS_H

#ifdef BOOST

#include <boost/math/constants/constants.hpp>

#endif

template<typename Scalar>
class constants {
public:
#ifdef BOOST
    static constexpr Scalar PI = boost::math::constants::pi<Scalar>();
#else
    static constexpr Scalar PI = static_cast<Scalar>(3.1415926535897932384626l);
#endif
};

template<typename Scalar>
constexpr Scalar constants<Scalar>::PI;

#endif