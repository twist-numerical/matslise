#ifndef SCHRODINGER_LEGENDRE_H
#define SCHRODINGER_LEGENDRE_H

#include "../matslise.h"

namespace legendre {
    template<typename Scalar>
    struct LegendreData {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> weights;
        Eigen::Array<Scalar, Eigen::Dynamic, 1> nodes;

        static const LegendreData<Scalar> *initData();

        static const LegendreData<Scalar> *getData(int n) {
            static const LegendreData<Scalar> *data = initData();
            int k = (n + 1) / 2;
            return data + k;
        }
    };


    template<class D, class Scalar=double>
    D *getCoefficients(int n, const std::function<D(Scalar)> &V, const Scalar &a, const Scalar &b) {
        const LegendreData<Scalar> *data = LegendreData<Scalar>::getData(n);

        Scalar m = (a + b) / 2;
        Scalar h = (b - a) / 2;
        Eigen::Array<D, Eigen::Dynamic, 1> result = (m + data->nodes * h).unaryExpr(V);
        D *coeffs = new D[n];
        Scalar H(1);
        for (int i = 0; i < n; ++i) {
            coeffs[i] = data->weights(i, 0) * result[0];
            for (int j = 1; j < data->weights.rows(); ++j)
                coeffs[i] += data->weights(i, j) * result[j];
            coeffs[i] /= H;
            H *= h;
        }
        return coeffs;
    }
}


#endif //SCHRODINGER_LEGENDRE_H
