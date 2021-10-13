#define INSTANTIATE_MATSLISE(Scalar)\
template class matslise::AbstractMatslise<Scalar>;\
template class matslise::Matslise<Scalar>;\
template class matslise::MatsliseHalf<Scalar>;\
template class matslise::PeriodicMatslise<Scalar>;


#include "../util/instantiate.h"
