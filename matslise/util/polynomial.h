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

        explicit Polynomial(const Scalar &v = 0) {
            coefficients = decltype(coefficients)::Zero();
            coefficients[0] = v;
        }

        constexpr Polynomial(std::initializer_list<Scalar> l) {
            assert(l.size() <= degree+1);
            Eigen::Index i = 0;
            for (const Scalar &s: l)
                coefficients[i++] = s;
            while (i <= degree) {
                coefficients[i++] = 0;
            }
        }

        explicit Polynomial(const Eigen::Array<Scalar, degree_ + 1, 1> &v) {
            coefficients = v;
        }

        template<int otherDegree, typename=typename std::enable_if<otherDegree <= degree, void>::type>
        Polynomial(const Polynomial<Scalar, otherDegree> &poly) {
            coefficients.template topRows<otherDegree + 1>() = poly.coefficients;
            constexpr int bottom = degree - otherDegree;
            if constexpr(bottom > 0)
                coefficients.template bottomRows<bottom>() = Eigen::Array<Scalar, bottom, 1>::Zero();
        }

#define POLY_INPLACE_OPERATOR(op, tr, rhs, dataField, sass) \
        typename std::enable_if<sass, Polynomial<Scalar, degree>>::type &operator op##=(const rhs &other) { \
            coefficients.template topRows<tr + 1>() op##= dataField; \
            return *this; \
        }

#define POLY_BINARY_OPERATOR(op, rhs) \
        Polynomial<Scalar, degree> operator op(const rhs &other) const { \
            Polynomial<Scalar, degree> r{*this}; \
            r op##= other; \
            return r; \
        }

#define POLY_BINARY_REV_OPERATOR(op, rhs) \
        friend Polynomial<Scalar, degree> operator op(const rhs &left, const Polynomial<Scalar, degree> &right) { \
            return right op left; \
        }

#define POLY_SCALAR_OPERATOR(op, tr) \
        POLY_INPLACE_OPERATOR(op, tr, Scalar, other, true); \
        POLY_BINARY_OPERATOR(op, Scalar); \

        template<int otherDegree>
        using OtherPolynomial = Polynomial<Scalar, otherDegree>;

#define POLY_SUM_OPERATOR(op) \
        template<int od> \
        POLY_INPLACE_OPERATOR(op, od, OtherPolynomial<od>, other.coefficients, od <= degree); \
        template<int od> \
        POLY_BINARY_OPERATOR(op, OtherPolynomial<od>);

        POLY_SUM_OPERATOR(+)

        template<int od>
        POLY_BINARY_REV_OPERATOR(+, OtherPolynomial<od>);

        POLY_SUM_OPERATOR(-)


        template<int od>
        friend Polynomial<Scalar, degree>
        operator-(const OtherPolynomial<od> &left, const Polynomial<Scalar, degree> &right) {
            return -right + left;
        }

        POLY_SCALAR_OPERATOR(+, 0)

        POLY_BINARY_REV_OPERATOR(+, Scalar);

        POLY_SCALAR_OPERATOR(-, 0)

        friend Polynomial<Scalar, degree> operator-(const Scalar &left, const Polynomial<Scalar, degree> &right) {
            return -right + left;
        }

        POLY_SCALAR_OPERATOR(*, degree)

        POLY_BINARY_REV_OPERATOR(*, Scalar);

        POLY_SCALAR_OPERATOR(/, degree)

        Polynomial<Scalar, degree> operator-() const {
            return Polynomial<Scalar, degree>{-coefficients};
        }

        Polynomial<Scalar, degree> operator+() const {
            return Polynomial<Scalar, degree>{+coefficients};
        }

        template<int otherDegree>
        Polynomial<Scalar, degree + otherDegree>
        operator*(const Polynomial<Scalar, otherDegree> &other) const {
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

        Scalar integral(Scalar x) const {
            Scalar r = 0;
            for (int i = degree; i >= 0; --i) {
                Scalar f = coefficients[i] / (i + 1);
                r = r * x + f;
            }
            return r * x;
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
