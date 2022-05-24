#define INSTANTIATE_MATSLISE(Scalar)\
template class matslise::LiouvilleTransformation<Scalar>;\
template class matslise::SturmLiouville<Scalar>;

#include "../util/instantiate.h"
