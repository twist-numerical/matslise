#ifndef MATSLISE_EIGEN_H
#define MATSLISE_EIGEN_H


#if __GNUG__ >= 8
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#pragma GCC diagnostic pop

#elif _MSC_VER
#pragma warning(push, 0)
MATSLISE_EIGEN_INCLUDE
#pragma warning(pop)

#else
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#endif

#endif //MATSLISE_EIGEN_H
