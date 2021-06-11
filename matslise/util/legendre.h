#ifndef SCHRODINGER_LEGENDRE_H
#define SCHRODINGER_LEGENDRE_H

#include "./eigen.h"

namespace matslise::legendre {
    template<typename Scalar, int n>
    const Eigen::Array<Scalar, n, n> weights;
    template<typename Scalar, int n>
    const Eigen::Array<Scalar, n, 1> nodes;

    template<int n, class D, class Scalar=double>
    std::array<D, n> getCoefficients(const std::function<D(Scalar)> &V, const Scalar &a, const Scalar &b) {
        constexpr const int k = n + n % 2;
        Scalar m = (a + b) / 2;
        Scalar h = (b - a) / 2;
        Eigen::Array<D, k, 1> result = (m + nodes<Scalar, k> * h).unaryExpr(V);
        std::array<D, n> coeffs;
        Scalar H(1);
        for (int i = 0; i < n; ++i) {
            coeffs[i] = weights<Scalar, k>(i, 0) * result[0];
            for (int j = 1; j < weights<Scalar, k>.rows(); ++j)
                coeffs[i] += weights<Scalar, k>(i, j) * result[j];
            coeffs[i] /= H;
            H *= h;
        }
        return coeffs;
    }
}

#include "./legendre_data.h"

#endif //SCHRODINGER_LEGENDRE_H
