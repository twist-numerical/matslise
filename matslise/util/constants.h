#ifndef MATSLISE_CONSTANTS_H
#define MATSLISE_CONSTANTS_H

#include <math.h>

template<typename Scalar>
class constants {
public:
    static const Scalar PI;
    static const Scalar SQRT1_2;
};

template<typename Scalar>
inline const Scalar constants<Scalar>::SQRT1_2 = sqrt(Scalar(.5));

#ifdef BOOST
#include <boost/math/constants/constants.hpp>
template<typename Scalar>
inline const Scalar constants<Scalar>::PI = boost::math::constants::pi<Scalar>();
#else
template<typename Scalar>
inline const Scalar constants<Scalar>::PI = static_cast<Scalar>(3.1415926535897932384626l);
#endif


#endif