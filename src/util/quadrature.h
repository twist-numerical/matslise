#ifndef MATSLISE_QUADRATURE_H
#define MATSLISE_QUADRATURE_H

#include <stack>
#include <tuple>
#include "eigen.h"


namespace quadrature {
    namespace lobatto {
        template<typename Scalar>
        Eigen::Array<Scalar, Eigen::Dynamic, 1> grid(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) {
            Eigen::Index n = x.size();
            Eigen::Array<Scalar, Eigen::Dynamic, 1> fine(3 * n - 2);
            Scalar sqrt5 = sqrt(static_cast<Scalar>(5));
            for (Eigen::Index i = 0; i < n - 1; ++i) {
                Scalar a = x[i], b = x[i + 1];
                Scalar m = (a + b) / 2, h = (b - a) / 2;
                fine[3 * i] = a;
                fine[3 * i + 1] = m - h / sqrt5;
                fine[3 * i + 2] = m + h / sqrt5;
            }
            fine[3 * (n - 1)] = x[n - 1];

            return fine;
        }

        template<typename Scalar>
        Scalar quadrature(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x,
                          const Eigen::Array<Scalar, Eigen::Dynamic, 1> &f) {
            Eigen::Index n = x.size();
            Scalar result = 0;

            for (Eigen::Index i = 0; i < n - 1; i += 3)
                result += (x[i + 3] - x[i]) / 2 * ((f[i] + f[i + 3]) / 6 + (f[i + 1] + f[i + 2]) * 5 / 6);

            return result;
        }

        template<typename Scalar>
        Scalar quadrature(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x,
                          const Eigen::Array<Scalar, Eigen::Dynamic, 1> &y,
                          const Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> &f) {
            Eigen::Array<Scalar, Eigen::Dynamic, 1> result = Eigen::Array<Scalar, Eigen::Dynamic, 1>::Zero(x.size());

            for (Eigen::Index i = 0; i < y.size() - 1; i += 3)
                result +=
                        (y[i + 3] - y[i]) / 2 * ((f.col(i) + f.col(i + 3)) / 6 + (f.col(i + 1) + f.col(i + 2)) * 5 / 6);

            return lobatto::quadrature<Scalar>(x, result);
        }
    }

    namespace gauss_konrod {
        template<typename Scalar, typename Value=Scalar>
        inline std::pair<Value, Scalar> applyGaussKonrod(
                const std::function<Eigen::Array<Value, Eigen::Dynamic, 1>(
                        const Eigen::Array<Scalar, Eigen::Dynamic, 1> &)> &f, Scalar a, Scalar b,
                const std::function<Scalar(const Value &)> &error);

        template<typename Scalar, typename Value=Scalar>
        inline std::pair<Value, Scalar> applyGaussKonrod(
                const std::function<Value(const Scalar &)> &f, Scalar a, Scalar b,
                const std::function<Scalar(const Value &)> &error);

        template<typename Scalar, typename Value=Scalar, bool bulk = false>
        Value adaptive(
                const std::function<typename std::conditional<bulk, Eigen::Array<Value, Eigen::Dynamic, 1>, Value>::type(
                        const typename std::conditional<bulk, Eigen::Array<Scalar, Eigen::Dynamic, 1>, Scalar>::type &)> &f,
                Scalar a, Scalar b, const Scalar &tolerance,
                const std::function<Scalar(const Value &)> &error = [](const Value &x) -> Scalar {
                    return abs(x);
                });
    }
}

#endif