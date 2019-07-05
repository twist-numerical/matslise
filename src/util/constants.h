#ifndef MATSLISE_CONSTANTS_H
#define MATSLISE_CONSTANTS_H

#ifdef BOOST

#include <boost/math/constants/constants.hpp>

#endif

template<typename Scalar>
class constants {
public:
    inline static Scalar pi() {
#ifdef BOOST
        return boost::math::constants::pi<Scalar>();
#else
        return (Scalar) M_PI;
#endif
    }
};


#endif