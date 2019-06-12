//
// Created by toon on 6/6/19.
//

#ifndef MATSLISE_INSTANTIATE_H
#define MATSLISE_INSTANTIATE_H

#ifndef INSTANTIATE_MORE
#define INSTANTIATE_MORE(Scalar)
#endif

#define INSTANTIATE_MATSLISE(Scalar)\
template class matslise::Matslise<Scalar>;\
template class matslise::HalfRange<Scalar>;\
template class matslise::Matscs<Scalar>;\
template class matslise::SE2D<Scalar>;\
INSTANTIATE_MORE(Scalar)

INSTANTIATE_MATSLISE(double)

INSTANTIATE_MATSLISE(long double)

#ifdef MATSLISE_float128

#include <boost/multiprecision/float128.hpp>

INSTANTIATE_MATSLISE(boost::multiprecision::float128)

#endif

#endif //MATSLISE_INSTANTIATE_H
