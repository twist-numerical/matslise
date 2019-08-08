#ifndef MATSLISE_EIGEN_H
#define MATSLISE_EIGEN_H

#if __GNUG__ >= 8
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/KroneckerProduct>
#pragma GCC diagnostic pop
#else
#pragma warning(push, 0)
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/KroneckerProduct>
#pragma warning(pop)
#endif

#endif //MATSLISE_EIGEN_H
