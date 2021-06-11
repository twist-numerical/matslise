#define INSTANTIATE_MATSLISE(Scalar)\
template class matslise::Matslise<Scalar>;\
template class matslise::MatsliseHalf<Scalar>;


#include "../util/instantiate.h"
