//
// Created by toon on 5/19/20.
//

#ifndef MATSLISE_SECTOR_CORRECTION_POTENTIAL_H
#define MATSLISE_SECTOR_CORRECTION_POTENTIAL_H

#include "../util/calculateEta.h"

#define MATSLISE_2HMAX_delta (MATSLISE_HMAX_delta*2-1)

template<typename Scalar>
std::complex<Scalar> sinhc(const std::complex<Scalar> &x) {
    if (abs(x) < .5) {
        // https://en.wikipedia.org/wiki/Sinhc_function#Pad%C3%A9_approximation
        std::complex<Scalar> x2 = x * x;
        return (Scalar(1) + x2 * (Scalar(53272705) / Scalar(360869676) +
                                  x2 * (Scalar(38518909) / Scalar(7217393520) +
                                        x2 * (Scalar(269197963) / Scalar(3940696861920) +
                                              x2 * Scalar(4585922449) / Scalar(15605159573203200))))
               ) / (Scalar(1) + x2 * (Scalar(-2290747) / Scalar(120289892) +
                                      x2 * (Scalar(1281433) / Scalar(7217393520) +
                                            x2 * (Scalar(-560401) / Scalar(562956694560) +
                                                  x2 * Scalar(1029037) / Scalar(346781323848960)))));
    } else
        return sinh(x) / x;
}

template<typename Scalar>
Eigen::Array<Scalar, 3, MATSLISE_HMAX_delta>
eta_integrals(const Scalar &delta, const Scalar &dZ) {
    Eigen::Array<Scalar, 3, MATSLISE_HMAX_delta> I = Eigen::Array<Scalar, 3, MATSLISE_HMAX_delta>::Zero();

    Scalar *eta = calculateEta<Scalar>(delta * delta * dZ, 2);

    I(0, 0) = (delta * eta[0] * eta[1] + delta) / 2;
    I(1, 1) = (eta[0] * eta[0] - 1) / (2 * dZ);
    I(2, 2) = (delta * eta[0] * eta[1] - delta) / (2 * dZ);

    Scalar dn = 1;
    for (int n = 1; n < MATSLISE_HMAX_delta; ++n) {
        if (n > 2) {
            I(2, n) = dn / delta * I(2, 2) - (n - 2) / (2 * dZ) * (I(1, n - 2) - dn / (n - 1));
        }
        if (n > 1) {
            I(1, n) = dn * I(1, 1) - 1 / (2 * dZ) * ((n - 1) * I(0, n - 2) - dn);
        }
        dn *= delta;
        I(0, n) = dn * I(0, 0) - Scalar(n) / 2 * (I(1, n) + dn * delta / (n + 1));
    }

    std::cout << "\n****\n" << delta << ", " << dZ << std::endl;
    std::cout << "\n" << I << std::endl;

    return I;
}

template<typename Scalar>
Eigen::Array<Eigen::Array<Scalar, MATSLISE_2HMAX_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta>
eta_integrals(const Scalar &delta, const Scalar &dZ1, const Scalar &dZ2) {
    Eigen::Array<Eigen::Array<Scalar, MATSLISE_2HMAX_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> I;
    for (int i = 0; i < MATSLISE_ETA_delta; ++i)
        for (int j = 0; j < MATSLISE_ETA_delta; ++j)
            I(i, j) = Eigen::Array<Scalar, MATSLISE_2HMAX_delta, 1>::Zero();

    Scalar Z1 = delta * delta * dZ1;
    Scalar Z2 = delta * delta * dZ2;
    Scalar *eta1 = calculateEta<Scalar>(Z1, 2);
    Scalar *eta2 = calculateEta<Scalar>(Z2, 2);
/*
    std::cout << "\n****\n" << delta << ", " << dZ1 << ", " << dZ2 << std::endl;
    std::cout << "\n****\n" << delta << ", " << dZ1 << ", " << dZ2 << std::endl;

    {
        Scalar denominator = 1 / (dZ1 - dZ2);
        I(0, 0)(0) = denominator * delta * (dZ1 * eta1[1] * eta2[0] - dZ2 * eta1[0] * eta2[1]);
        I(1, 0)(1) = denominator * (eta1[0] * eta2[0] - Z2 * eta1[1] * eta2[1] - 1);
        I(0, 1)(1) = -denominator * (eta1[0] * eta2[0] - Z1 * eta1[1] * eta2[1] - 1);
        I(1, 1)(2) = denominator * delta * (eta1[0] * eta2[1] - eta1[1] * eta2[0]);

        Scalar dn = 1;
        for (int n = 1; n < MATSLISE_2HMAX_delta; ++n) {
            if (n > 2) {
                I(1, 1)(n) = dn / delta * I(1, 1)(2) - denominator * (n - 2) * (I(0, 1)(n - 2) - I(1, 0)(n - 2));
            }
            if (n > 1) {
                I(1, 0)(n) =
                        dn * I(1, 0)(1) - denominator * (n - 1) * (I(0, 0)(n - 2) - dZ2 * I(1, 1)(n) - dn / (n - 1));
                I(0, 1)(n) =
                        dn * I(0, 1)(1) + denominator * (n - 1) * (I(0, 0)(n - 2) - dZ1 * I(1, 1)(n) - dn / (n - 1));
            }
            dn *= delta;
            I(0, 0)(n) = dn * I(0, 0)(0) - denominator * n * (dZ1 * I(1, 0)(n) - dZ2 * I(0, 1)(n));
        }

        Eigen::Array<Scalar, 4, MATSLISE_2HMAX_delta> coutI;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                coutI.row(i * 2 + j) = I(i, j);
        std::cout << coutI << "\n" << std::endl;

        for (int i = 0; i < MATSLISE_ETA_delta; ++i)
            for (int j = 0; j < MATSLISE_ETA_delta; ++j)
                I(i, j) = Eigen::Array<Scalar, MATSLISE_2HMAX_delta, 1>::Zero();
    }*/


    {
        Scalar denominator = 1 / (dZ1 - dZ2);
        Scalar &I000 = I(0, 0)(0) = denominator * delta * (dZ1 * eta1[1] * eta2[0] - dZ2 * eta1[0] * eta2[1]);
        Scalar &I101 = I(1, 0)(1) = denominator * (eta1[0] * eta2[0] - Z2 * eta1[1] * eta2[1] - 1);
        Scalar &I011 = I(0, 1)(1) = -denominator * (eta1[0] * eta2[0] - Z1 * eta1[1] * eta2[1] - 1);
        Scalar &I112 = I(1, 1)(2) = denominator * delta * (eta1[0] * eta2[1] - eta1[1] * eta2[0]);

        int lastN = MATSLISE_2HMAX_delta - 3;
        Scalar dn = std::pow(delta, lastN);
        for (int n = lastN; n > 1; --n) {
            Scalar x1 = (I101 * dn - I(1, 0)(n + 1));
            Scalar x2 = (I011 * dn - I(0, 1)(n + 1));
            I(0, 0)(n - 1) = (x1 * dZ1 + x2 * dZ2 + dn) / n;
            I(1, 1)(n + 1) = (x1 + x2) / n;

            Scalar y1 = I000 * dn - I(0, 0)(n);
            Scalar y2 = I112 * dn - I(1, 1)(n + 2);
            I(0, 1)(n) = (y1 + y2 * dZ1) / n;
            I(1, 0)(n) = (y1 + y2 * dZ2) / n;

            dn /= delta;
        }

        /*
        Eigen::Array<Scalar, 4, MATSLISE_2HMAX_delta> coutI;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                coutI.row(i * 2 + j) = I(i, j);
        std::cout << coutI << "\n" << std::endl;
         */
    }


    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        if (i > 1) {
            for (int j = 0; j < 2; ++j)
                I(i, j).template bottomRows<MATSLISE_2HMAX_delta - 2>()
                        = ((I(i - 2, j) - (2 * i - 3) * I(i - 1, j)) / dZ1)
                        .template topRows<MATSLISE_2HMAX_delta - 2>();
        }
        for (int j = 2; j < MATSLISE_ETA_delta; ++j) {
            I(i, j).template bottomRows<MATSLISE_2HMAX_delta - 2>()
                    = ((I(i, j - 2) - (2 * j - 3) * I(i, j - 1)) / dZ2)
                    .template topRows<MATSLISE_2HMAX_delta - 2>();
        }
    }

    /*{
        Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, MATSLISE_2HMAX_delta> coutI;
        for (int i = 0; i < MATSLISE_ETA_delta; ++i)
            for (int j = 0; j < MATSLISE_ETA_delta; ++j)
                coutI.row(i * MATSLISE_ETA_delta + j) = I(i, j);
        std::cout << coutI << std::endl;
    }*/


    return I;
}

template<typename Scalar>
Eigen::Array<std::complex<Scalar>, MATSLISE_HMAX_delta * 2 - 1, 1>
exp_integrals(const Scalar &delta, const std::complex<Scalar> &theta) {
    Eigen::Array<std::complex<Scalar>, MATSLISE_HMAX_delta * 2 - 1, 1> integrals
            = Eigen::Array<Scalar, MATSLISE_HMAX_delta * 2 - 1, 1>::Zero();

    // std::cout << "θδ(" << theta << ", " << delta << ")" << std::endl;
    if (abs(theta) < Scalar(5) || true) {
        Scalar d = std::pow(delta, MATSLISE_HMAX_delta * 2 - 2);
        std::complex<Scalar> etd = exp(theta * delta);
        for (int i = MATSLISE_HMAX_delta * 2 - 2; i > 0; --i) {
            // d = delta**i
            integrals(i - 1) = (d * etd - theta * integrals(i)) / Scalar(i);
            d /= delta;
        }
    } else {
        integrals(0) = (exp(theta * delta) - Scalar(1)) / theta;
        Scalar d = 1;
        for (int i = 1; i < MATSLISE_HMAX_delta - 1; ++i) {
            d *= delta;
            integrals(i) = d / theta * exp(theta * delta) - Scalar(i) / theta * integrals(i - 1);
        }
    }
    return integrals;
}

/*
template<typename Scalar>
Eigen::Array<Scalar, MATSLISE_N, 1> vbar_formulas(
        const Eigen::Array<Scalar, MATSLISE_HMAX_delta, MATSLISE_ETA_delta> &y, const Scalar &delta, const Scalar &dZ) {

    Eigen::Array<Eigen::Array<Scalar, MATSLISE_ETA_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> eta = eta_integrals(delta, dZ);


    Eigen::Array<Scalar, MATSLISE_N, 1> quadratures = Eigen::Array<Scalar, MATSLISE_N, 1>::Zero();


    std::cout << "\n --- quad " << std::endl;
    std::cout << "\n" << quadratures.transpose() << std::endl;

    static Eigen::Matrix<Scalar, MATSLISE_N, MATSLISE_N> legendrePolynomials = (
            Eigen::Matrix<Scalar, MATSLISE_N, MATSLISE_N>()
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

    Eigen::Matrix<Scalar, MATSLISE_N, 1> legendreScaling;
    Scalar d = 1;
    for (int i = 0; i < MATSLISE_N; ++i, d /= delta)
        legendreScaling(i) = d;
    return (legendrePolynomials * legendreScaling.asDiagonal() * quadratures.matrix()).array();
}
*/

template<typename Scalar>
Eigen::Array<Scalar, MATSLISE_N, 1> vbar_formulas(
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &y1,
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &y2,
        const Scalar &delta, const Scalar &dZ1, const Scalar &dZ2) {

    Eigen::Array<Eigen::Array<Scalar, MATSLISE_2HMAX_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta>
            eta = eta_integrals<Scalar>(delta, dZ1, dZ2);

    Eigen::Array<Scalar, MATSLISE_N, 1> quadratures = Eigen::Array<Scalar, MATSLISE_N, 1>::Zero();
/*
    std::cout << "\ny1: " << std::endl;
    std::cout << y1 << std::endl;
    std::cout << "\ny2: " << std::endl;
    std::cout << y2 << std::endl;
*/
    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        for (int j = 0; j < MATSLISE_ETA_delta; ++j) {
            for (int ni = 0; ni < MATSLISE_HMAX_delta; ++ni) {
                for (int nj = 0; nj < MATSLISE_HMAX_delta; ++nj) {
                    for (int k = 0; k < MATSLISE_N && ni + nj + k < MATSLISE_2HMAX_delta; ++k) {
                        quadratures(k) += eta(i, j)(ni + nj + k) * y1(i, ni) * y2(j, nj);
                        /*if (eta(i, j)(ni + nj + k) != 0 || y1(i, ni) != 0 || y2(j, nj) != 0) {
                            if (eta(i, j)(ni + nj + k) * y1(i, ni) * y2(j, nj) == 0) {
                                std::cout << eta(i, j)(ni + nj + k) << ", "
                                          << y1(i, ni) << ", " << y2(j, nj) << std::endl;
                                std::cout << "error!" << std::endl;
                            }
                        }*/
                    }
                }
            }
        }
    }
/*
    std::cout << "\n --- quad " << std::endl;
    std::cout << "\n" << quadratures.transpose() << std::endl;
*/
    static Eigen::Matrix<Scalar, MATSLISE_N, MATSLISE_N> legendrePolynomials = (
            Eigen::Matrix<Scalar, MATSLISE_N, MATSLISE_N>()
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

    Eigen::Matrix<Scalar, MATSLISE_N, 1> legendreScaling;
    Scalar d = 1;
    for (int i = 0; i < MATSLISE_N; ++i, d /= delta)
        legendreScaling(i) = d;
    return (legendrePolynomials * legendreScaling.asDiagonal() * quadratures.matrix()).array();
}

#endif //MATSLISE_SECTOR_CORRECTION_POTENTIAL_H