#ifndef INSTANTIATE_MORE
#define INSTANTIATE_MORE(Scalar)
#endif

#ifndef INSTANTIATE_MATSLISE
#define INSTANTIATE_MATSLISE(Scalar)
#endif

/*
#ifndef INSTANTIATE_MATSLISE
#define INSTANTIATE_MATSLISE(Scalar)\
template class matslise::Matslise<Scalar>;\
template class matslise::MatsliseHalf<Scalar>;\
template class matslise::Matscs<Scalar>;\
template class matslise::Matslise2DSector<Scalar>;\
template class matslise::Matslise2D<Scalar>;\
template class matslise::Matslise2DHalf<Scalar>;\
template class matslise::Matslise3DSector<Scalar>;\
template class matslise::Matslise3D<Scalar>;\
template class matslise::MatsliseNDSector<Scalar>;\
template class matslise::MatsliseND<Scalar, matslise::Matslise2DSector<Scalar>>;\
template class matslise::MatsliseND<Scalar, matslise::Matslise3DSector<Scalar>>;\
INSTANTIATE_MORE(Scalar)
#endif
 */

#define INSTANTIATE_ALL(Scalar) \
INSTANTIATE_MATSLISE(Scalar) \
INSTANTIATE_MORE(Scalar)

INSTANTIATE_ALL(double)

#ifdef MATSLISE_LONG_DOUBLE

INSTANTIATE_ALL(long double)

#endif

#ifdef MATSLISE_QUADMATH

#include <boost/multiprecision/float128.hpp>

INSTANTIATE_ALL(boost::multiprecision::float128)

#endif
