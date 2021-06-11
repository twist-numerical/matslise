#define INSTANTIATE_MATSLISE(Scalar)\
template class matslise::Matslise2DSector<Scalar>;\
template class matslise::Matslise2D<Scalar>;\
template class matslise::Matslise2DHalf<Scalar>;\
template class matslise::MatsliseND<Scalar, matslise::Matslise2DSector<Scalar>>;

#include "../util/instantiate.h"
