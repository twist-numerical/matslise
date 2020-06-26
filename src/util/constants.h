#ifndef MATSLISE_CONSTANTS_H
#define MATSLISE_CONSTANTS_H

template<typename Scalar>
class constants {
public:
    static const Scalar PI;
};

#ifdef BOOST
#include <boost/math/constants/constants.hpp>

template<typename Scalar>
inline const Scalar constants<Scalar>::PI = boost::math::constants::pi<Scalar>();
#else
inline const Scalar constants<Scalar>::PI  = static_cast<Scalar>(3.1415926535897932384626l);
#endif


#endif