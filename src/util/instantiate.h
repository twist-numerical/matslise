//
// Created by toon on 6/6/19.
//

#ifndef MATSLISE_INSTANTIATE_H
#define MATSLISE_INSTANTIATE_H

#define INSTANTIATE_MATSLISE(Scalar) \
template class matslise::Matslise<Scalar>; \
template class matslise::HalfRange<Scalar>; \
template class matslise::Matscs<Scalar>;

INSTANTIATE_MATSLISE(double)

INSTANTIATE_MATSLISE(long double)

#ifdef BOOST

#include <boost/multiprecision/float128.hpp>

INSTANTIATE_MATSLISE(boost::multiprecision::float128)

#endif

#endif //MATSLISE_INSTANTIATE_H
