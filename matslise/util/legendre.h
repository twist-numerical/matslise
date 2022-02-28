#ifndef SCHRODINGER_LEGENDRE_H
#define SCHRODINGER_LEGENDRE_H

#include "./eigen.h"

namespace matslise::legendre {
    template<typename Scalar, int n>
    const Eigen::Array<Scalar, n, n> weights;
    template<typename Scalar, int n>
    const Eigen::Array<Scalar, n, 1> nodes;

    template<typename Scalar>
    inline const Eigen::Matrix<Scalar, 16, 16> polynomials = (
            Eigen::Matrix<Scalar, 16, 16>()
                    << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    -1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, -6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    -1, 12, -30, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, -20, 90, -140, 70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    -1, 30, -210, 560, -630, 252, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, -42, 420, -1680, 3150, -2772, 924, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    -1, 56, -756, 4200, -11550, 16632, -12012, 3432, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, -72, 1260, -9240, 34650, -72072, 84084, -51480, 12870, 0, 0, 0, 0, 0, 0, 0,
                    -1, 90, -1980, 18480, -90090, 252252, -420420, 411840, -218790, 48620, 0, 0, 0, 0, 0, 0,
                    1, -110, 2970, -34320, 210210, -756756, 1681680, -2333760, 1969110, -923780, 184756, 0, 0, 0, 0, 0,
                    -1, 132, -4290, 60060, -450450, 2018016, -5717712, 10501920, -12471030, 9237800, -3879876, 705432, 0, 0, 0, 0,
                    1, -156, 6006, -100100, 900900, -4900896, 17153136, -39907296, 62355150, -64664600, 42678636, -16224936, 2704156, 0, 0, 0,
                    -1, 182, -8190, 160160, -1701700, 11027016, -46558512, 133024320, -261891630, 355655300, -327202876, 194699232, -67603900, 10400600, 0, 0,
                    1, -210, 10920, -247520, 3063060, -23279256, 116396280, -399072960, 960269310, -1636014380, 1963217256, -1622493600, 878850700, -280816200, 40116600, 0,
                    -1, 240, -14280, 371280, -5290740, 46558512, -271591320, 1097450640, -3155170590, 6544057520, -9816086280, 10546208400, -7909656300, 3931426800, -1163381400, 155117520
    ).finished();

    template<int n, typename Scalar, typename D=Scalar>
    class Legendre {
        static const constexpr int k = n + n % 2;
        const Scalar h;
        const Eigen::Array<D, k, 1> fs;
    public:
        template<typename F>
        Legendre(const F &f, const Scalar &a, const Scalar &b)
                : h((b - a) / 2), fs(((a + b) / 2 + nodes<Scalar, k> * h).unaryExpr(f)) {}

        D integrate(int i = 0) {
            D result = weights<Scalar, k>(i, 0) * fs[0];
            for (int j = 1; j < weights<Scalar, k>.rows(); ++j)
                result += weights<Scalar, k>(i, j) * fs[j];
            return result;
        }

        std::array<D, n> getCoefficients() {
            std::array<D, n> coeffs;
            Scalar H(1);
            for (int i = 0; i < n; ++i) {
                coeffs[i] = integrate(i) / H;
                H *= h;
            }
            return coeffs;
        }

        std::array<D, n> asPolynomial() {
            static_assert(n <= 16);

            std::array<D, n> coeffs = getCoefficients();
            std::array<D, n> polynomial{};
            Scalar H = 1;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j <= i; ++j) {
                    polynomial[j] += (H * polynomials<Scalar>(i, j)) * coeffs[i];
                }
                H *= 2 * h;
            }
            return polynomial;
        }
    };

    template<int n, class Scalar, class D = Scalar, typename F>
    std::array<D, n> getCoefficients(const F &f, const Scalar &a, const Scalar &b) {
        return Legendre<n, Scalar, D>{f, a, b}.getCoefficients();
    }
}

#include "./legendre_data.h"

#endif //SCHRODINGER_LEGENDRE_H
