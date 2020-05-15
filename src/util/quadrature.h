#ifndef MATSLISE_QUADRATURE_H
#define MATSLISE_QUADRATURE_H

#include <stack>
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

    namespace trapezoidal {
        template<typename Scalar, typename Value=Scalar>
        Value adaptive(const std::function<Value(Scalar)> &f, Scalar a, Scalar b,
                       const Scalar &tolerance,
                       const std::function<Scalar(const Value &)> error = [](const Value &x) -> Scalar {
                           return abs(x);
                       }) {
            std::stack<std::tuple<Scalar, Value, Scalar, Value>> todo;
            Scalar l = b - a;
            todo.push({a, f(a), b, f(b)});
            Value fa, fb;
            Value result;
            bool resultInitialized = false;
            while (!todo.empty()) {
                std::tie(a, fa, b, fb) = todo.top();
                todo.pop();
                Scalar c = (a + b) / 2;
                Value fc = f(c);
                Value estimate = (b - a) / 6 * (fa + 4 * fc + fb);
                if (error(estimate - (b - a) / 4 * (fa + 2 * fc + fb)) < tolerance / l * (b - a)) {
                    if (resultInitialized)
                        result += estimate;
                    else {
                        result = estimate;
                        resultInitialized = true;
                    }
                } else {
                    todo.push({a, fa, c, fc});
                    todo.push({c, fc, b, fb});
                }
            }
            return result;
        };
    }
}

#endif