#ifndef MATSLISE_INSTANTIATE_H
#define MATSLISE_INSTANTIATE_H

#ifndef INSTANTIATE_MORE
#define INSTANTIATE_MORE(Scalar)
#endif

#ifndef INSTANTIATE_MATSLISE
#define INSTANTIATE_MATSLISE(Scalar)\
template class matslise::Matslise<Scalar>;\
template class matslise::MatsliseHalf<Scalar>;\
template class matslise::Matscs<Scalar>;\
template class matslise::Matslise2D<Scalar>;\
template class matslise::Matslise2DHalf<Scalar>;\
template class matslise::Matslise3D<Scalar>;\
INSTANTIATE_MORE(Scalar)

#endif

INSTANTIATE_MATSLISE(double)

#ifdef MATSLISE_long_double

INSTANTIATE_MATSLISE(long double)

#endif

#ifdef MATSLISE_float128

#include <boost/multiprecision/float128.hpp>

INSTANTIATE_MATSLISE(boost::multiprecision::float128)

#endif

#endif //MATSLISE_INSTANTIATE_H
