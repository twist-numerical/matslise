#ifndef MATSLISE_POLYNOMIAL_H
#define MATSLISE_POLYNOMIAL_H

#include "./eigen.h"

namespace matslise {
    template<typename Scalar, int degree_>
    class Polynomial {
    public:
        template<typename, int> friend
        class Polynomial;

        static const constexpr int degree = degree_;
    private:
        Eigen::Array<Scalar, degree_ + 1, 1> coefficients;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Polynomial(const Scalar &v = 0) {
            coefficients = decltype(coefficients)::Zero();
            coefficients[0] = v;
        }

        Polynomial(const std::array<Scalar, degree + 1> &v) {
            for (unsigned int i = 0; i < v.size(); ++i)
                coefficients[i] = v[i];
        }

        template<int otherDegree>
        Polynomial(const Polynomial<Scalar, otherDegree> &poly) {
            static_assert(otherDegree <= degree);
            coefficients.template topRows<otherDegree + 1>() = poly.coefficients;
            constexpr int bottom = degree - otherDegree;
            if constexpr(bottom > 0)
                coefficients.template bottomRows<bottom>() = Eigen::Array<Scalar, bottom, 1>::Zero();
        }

        Polynomial<Scalar, degree> &operator*=(const Scalar &f) {
            coefficients *= f;
            return *this;
        }

#define POLY_SIMPLE_OPERATOR(op) \
    template<int otherDegree> \
    Polynomial<Scalar, degree> &operator op##=(const Polynomial<Scalar, otherDegree> &other) { \
        static_assert(otherDegree <= degree); \
        coefficients.template topRows<otherDegree + 1>() op##= other.coefficients; \
        return *this; \
    }

        POLY_SIMPLE_OPERATOR(+)

        POLY_SIMPLE_OPERATOR(-)

        template<int otherDegree>
        Polynomial<Scalar, degree + otherDegree> operator*(const Polynomial<Scalar, otherDegree> &other) const {
            Polynomial<Scalar, degree + otherDegree> r{};
            for (int i = 0; i <= degree; ++i)
                for (int j = 0; j <= otherDegree; ++j)
                    r.coefficients[i + j] += coefficients[i] * other.coefficients[j];
            return r;
        }

        bool operator==(const Polynomial<Scalar, degree> &other) const {
            return (coefficients == other.coefficients).all();
        }

        bool operator!=(const Polynomial<Scalar, degree> &other) const {
            return !operator==(other);
        }

        Scalar &operator[](int x) {
            return coefficients[x];
        }

        const Scalar &operator[](int x) const {
            return coefficients[x];
        }

        Scalar operator()(const Scalar &x) const {
            Scalar r = coefficients[degree];
            for (int i = degree - 1; i >= 0; --i)
                r = r * x + coefficients[i];
            return r;
        }

        template<int d = 1>
        Polynomial<Scalar, degree - d> derivative() const {
            static_assert(d >= 0);
            Polynomial<Scalar, degree - d> r{};

            for (int i = d; i <= degree; ++i) {
                Scalar f = coefficients[i];
                for (int j = i; j > i - d; --j) f *= j;
                r[i - d] = f;
            }
            return r;
        }

        template<int d = 1>
        Scalar derivative(Scalar x) const {
            static_assert(d >= 0);

            Scalar r = 0;
            for (int i = degree; i >= d; --i) {
                Scalar f = coefficients[i];
                for (int j = i; j > i - d; --j) f *= j;
                r = r * x + f;
            }
            return r;
        }

        Polynomial<Scalar, degree + 1> integral() const {
            Polynomial<Scalar, degree + 1> r{};
            for (int i = 0; i <= degree; ++i) {
                r[i + 1] = coefficients[i] / (i + 1);
            }
            return r;
        }

        Eigen::Index size() const {
            return coefficients.size();
        }
    };

#define POLY_FRIEND_OPERATOR(op) \
    template<typename Scalar, int lhsDegree, int rhsDegree> \
    Polynomial<Scalar, std::max(lhsDegree, rhsDegree)> operator op( \
            const Polynomial<Scalar, lhsDegree> &lhs, const Polynomial<Scalar, rhsDegree> &rhs) { \
        Polynomial<Scalar, std::max(lhsDegree, rhsDegree)> result{lhs}; \
        return result op##= rhs; \
    }

    POLY_FRIEND_OPERATOR(+)

    POLY_FRIEND_OPERATOR(-)
}

#endif //MATSLISE_POLYNOMIAL_H
