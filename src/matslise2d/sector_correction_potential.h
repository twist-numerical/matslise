//
// Created by toon on 5/19/20.
//

#ifndef MATSLISE_SECTOR_CORRECTION_POTENTIAL_H
#define MATSLISE_SECTOR_CORRECTION_POTENTIAL_H

#include "../util/calculateEta.h"

#define MATSLISE_INTEGRATE_delta (2*MATSLISE_HMAX_delta)


template<typename Scalar, Eigen::Index N>
Eigen::Array<Eigen::Array<Scalar, N, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> etaProduct(
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &u,
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &v) {

    Eigen::Array<Eigen::Array<Scalar, N, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> uv;

    for (int i = 0; i < MATSLISE_ETA_delta; ++i)
        for (int j = 0; j < MATSLISE_ETA_delta; ++j) {
            uv(i, j) = Eigen::Array<Scalar, N, 1>::Zero();

            for (int ni = 0; ni < MATSLISE_HMAX_delta; ++ni)
                for (int nj = 0; nj < MATSLISE_HMAX_delta && ni + nj < N; ++nj)
                    uv(i, j)(ni + nj) += u(i, ni) * v(j, nj);
        }
    return uv;
}

template<typename Scalar>
Scalar powInt(Scalar x, unsigned int n) {
    if (n == 0)
        return 1;
    Scalar y = 1;
    for (; n > 1; n >>= 1) {
        if (n & 1)
            y *= x;
        x *= x;
    }
    return x * y;
}


template<typename Scalar>
Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> integrateTaylorTaylor(
        const Scalar &delta,
        const Scalar &dZ1, const Scalar *eta1,
        const Scalar &dZ2, const Scalar *eta2,
        int degree
) {
    Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> result
            = Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>::Zero();
    Scalar x1n = 1;
    for (int i = 0; i <= degree; ++i) {
        Scalar x2n = 1;
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k < MATSLISE_INTEGRATE_delta; ++k)
                result(k) += eta1[i] * x1n * eta2[j] * x2n * powInt(delta, 2 * i + 2 * j + k + 1) /
                             (2 * i + 2 * j + k + 1);
            x2n *= dZ2;
        }
        x1n *= dZ1;
    }
    return result;
}

template<typename Scalar, int degree>
Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 2> integrateEtaTaylor(
        const Scalar &delta,
        const Scalar &dZ1, const Scalar &dZ2, const Scalar *taylor) {
    Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 2> I
            = Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 2>::Zero();

    Scalar *eta = calculateEta<Scalar>(delta * delta * dZ1, 2);

    Scalar xi_i = 0;
    Scalar eta0_i = 0;

    Scalar di = powInt(delta, MATSLISE_INTEGRATE_delta + 2 * degree - 1);
    for (int i = MATSLISE_INTEGRATE_delta + 2 * degree - 2; i > 0; --i) {
        Scalar prev_eta0_i = eta0_i;
        eta0_i = (di * eta[1] - xi_i) / i;
        //result(i, 1) = (di * eta[1] - result(i, 0)) / i;
        di /= delta;
        xi_i = (di * eta[0] - dZ1 * prev_eta0_i) / i;
        //result(i - 1, 0) = (di * eta[0] - dZ1 * result(i + 1, 1)) / i;
        for (int k = std::max(0, i - MATSLISE_INTEGRATE_delta + 1); k <= degree && 2 * k <= i; ++k) {
            Scalar dZ2k = powInt(dZ2, k);
            if (2 * k < i)
                I(i - 2 * k - 1, 0) += dZ2k * taylor[k] * xi_i;
            I(i - 2 * k, 1) += dZ2k * taylor[k] * eta0_i;
        }
    }
/*
    for (int k = 0; k <= degree; ++k) {
        std::cout << taylor[k] << ", ";
    }
    std::cout << std::endl;

    std::cout << "\n" << I.transpose() << std::endl;
*/
    return I;
}

template<typename Scalar, bool equal>
Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta>
eta_integrals(const Scalar &delta, const Scalar &dZ1, const Scalar &dZ2) {
    Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> I;
    for (int i = 0; i < MATSLISE_ETA_delta; ++i)
        for (int j = 0; j < MATSLISE_ETA_delta; ++j)
            I(i, j) = Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>::Zero();

    Scalar Z1 = delta * delta * dZ1;
    Scalar Z2 = delta * delta * dZ2;

    if (abs(Z1) < 0.5 && abs(Z2) < 0.5) {
        Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, 2, 2> step;
        step << integrateTaylorTaylor(delta, dZ1, CalculateEtaData<Scalar>::taylor_eta8, dZ2,
                                      CalculateEtaData<Scalar>::taylor_eta8, CalculateEtaData<Scalar>::taylor_degree),
                integrateTaylorTaylor(delta, dZ1, CalculateEtaData<Scalar>::taylor_eta8, dZ2,
                                      CalculateEtaData<Scalar>::taylor_eta9, CalculateEtaData<Scalar>::taylor_degree),
                integrateTaylorTaylor(delta, dZ1, CalculateEtaData<Scalar>::taylor_eta9, dZ2,
                                      CalculateEtaData<Scalar>::taylor_eta8, CalculateEtaData<Scalar>::taylor_degree),
                integrateTaylorTaylor(delta, dZ1, CalculateEtaData<Scalar>::taylor_eta9, dZ2,
                                      CalculateEtaData<Scalar>::taylor_eta9, CalculateEtaData<Scalar>::taylor_degree);

        /*
        std::cout << step(0, 0).transpose() << std::endl;
        std::cout << step(0, 1).transpose() << std::endl;
        std::cout << step(1, 0).transpose() << std::endl;
        std::cout << step(1, 1).transpose() << std::endl;
        */

        for (int i = 8; i >= 0; --i) {
            {
                for (int j = 0; j < 2; ++j) {
                    Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = step(0, j);
                    step(0, j)(MATSLISE_INTEGRATE_delta - 1) = 0;
                    step(0, j)(MATSLISE_INTEGRATE_delta - 2) = 0;
                    step(0, j).template topRows<MATSLISE_INTEGRATE_delta - 2>() =
                            dZ1 * step(1, j).template bottomRows<MATSLISE_INTEGRATE_delta - 2>()
                            + (2 * i + 1) * curr.template topRows<MATSLISE_INTEGRATE_delta - 2>();
                    step(1, j) = curr;
                }
            }
            Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, 1, 2> stepUp = step.row(0);
            for (int j = 8; j >= 0; --j) {
                Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = stepUp(0);
                stepUp(0)(MATSLISE_INTEGRATE_delta - 1) = 0;
                stepUp(0)(MATSLISE_INTEGRATE_delta - 2) = 0;
                stepUp(0).template topRows<MATSLISE_INTEGRATE_delta - 2>() =
                        dZ2 * stepUp(1).template bottomRows<MATSLISE_INTEGRATE_delta - 2>()
                        + (2 * j + 1) * curr.template topRows<MATSLISE_INTEGRATE_delta - 2>();
                stepUp(1) = curr;
                if (i < MATSLISE_ETA_delta && j < MATSLISE_ETA_delta)
                    I(i, j) = stepUp(0);
            }
            /*
            I(i, j).template bottomRows<MATSLISE_INTEGRATE_delta - 2>()
                    = ((I(i, j - 2) - (2 * j - 3) * I(i, j - 1)) / dZ2)
                    .template topRows<MATSLISE_INTEGRATE_delta - 2>();
             */
        }
/*
        std::cout << "\n****\n" << delta << ", " << dZ1 << ", " << dZ2 << std::endl;

        Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, MATSLISE_INTEGRATE_delta> coutI;
        for (int i = 0; i < MATSLISE_ETA_delta; ++i)
            for (int j = 0; j < MATSLISE_ETA_delta; ++j)
                coutI.row(i * MATSLISE_ETA_delta + j) = I(i, j);
        std::cout << coutI << std::endl;
*/
        return I;
    } else if (abs(Z1) < 0.5) {
        Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, 2, 2> step;
        Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 2> e8
                = integrateEtaTaylor<Scalar, CalculateEtaData<Scalar>::taylor_degree>(
                        delta, dZ2, dZ1, CalculateEtaData<Scalar>::taylor_eta8);
        step(0, 0) = e8.col(0);
        step(0, 1) = e8.col(1);
        Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 2> e9
                = integrateEtaTaylor<Scalar, CalculateEtaData<Scalar>::taylor_degree>(
                        delta, dZ2, dZ1, CalculateEtaData<Scalar>::taylor_eta9);
        step(1, 0) = e9.col(0);
        step(1, 1) = e9.col(1);

        for (int i = 8; i >= 0; --i) {
            {
                for (int j = 0; j < 2; ++j) {
                    Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = step(0, j);
                    step(0, j)(MATSLISE_INTEGRATE_delta - 1) = 0;
                    step(0, j)(MATSLISE_INTEGRATE_delta - 2) = 0;
                    step(0, j).template topRows<MATSLISE_INTEGRATE_delta - 2>() =
                            dZ1 * step(1, j).template bottomRows<MATSLISE_INTEGRATE_delta - 2>()
                            + (2 * i + 1) * curr.template topRows<MATSLISE_INTEGRATE_delta - 2>();
                    step(1, j) = curr;
                }
            }

            Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, 1, 2> stepUp = step.row(0);
            for (int j = 0; j < MATSLISE_ETA_delta; ++j) {
                if (i < MATSLISE_ETA_delta)
                    I(i, j) = stepUp(0);
                Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = stepUp(1);
                stepUp(1)(0) = 0;
                stepUp(1)(1) = 0;
                stepUp(1).template tail<MATSLISE_INTEGRATE_delta - 2>()
                        = ((stepUp(0) - (2 * j + 1) * curr) / dZ2)
                        .template head<MATSLISE_INTEGRATE_delta - 2>();
                stepUp(0) = curr;
            }
        }

        /*{
            std::cout << "\n****\n" << delta << ", " << dZ1 << ", " << dZ2 << std::endl;

            Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, MATSLISE_INTEGRATE_delta> coutI;
            for (int i = 0; i < MATSLISE_ETA_delta; ++i)
                for (int j = 0; j < MATSLISE_ETA_delta; ++j)
                    coutI.row(i * MATSLISE_ETA_delta + j) = I(i, j);
            std::cout << coutI << std::endl;
            std::cout << std::endl;
        }*/
        return I;
    } else if (abs(Z2) < 0.5) {
        Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, 2, 2> step;
        Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 2> e8
                = integrateEtaTaylor<Scalar, CalculateEtaData<Scalar>::taylor_degree>(
                        delta, dZ1, dZ2, CalculateEtaData<Scalar>::taylor_eta8);
        step(0, 0) = e8.col(0);
        step(1, 0) = e8.col(1);
        Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 2> e9
                = integrateEtaTaylor<Scalar, CalculateEtaData<Scalar>::taylor_degree>(
                        delta, dZ1, dZ2, CalculateEtaData<Scalar>::taylor_eta9);
        step(0, 1) = e9.col(0);
        step(1, 1) = e9.col(0);


        for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
            Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, 1, 2> stepUp = step.row(0);
            for (int j = 8; j >= 0; --j) {
                Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = stepUp(0);
                stepUp(0)(MATSLISE_INTEGRATE_delta - 1) = 0;
                stepUp(0)(MATSLISE_INTEGRATE_delta - 2) = 0;
                stepUp(0).template topRows<MATSLISE_INTEGRATE_delta - 2>() =
                        dZ2 * stepUp(1).template bottomRows<MATSLISE_INTEGRATE_delta - 2>()
                        + (2 * j + 1) * curr.template topRows<MATSLISE_INTEGRATE_delta - 2>();
                stepUp(1) = curr;
                if (i < MATSLISE_ETA_delta && j < MATSLISE_ETA_delta)
                    I(i, j) = stepUp(0);
            }
            {
                for (int j = 0; j < 2; ++j) {
                    Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = step(1, j);
                    step(1, j)(0) = 0;
                    step(1, j)(1) = 0;
                    step(1, j).template tail<MATSLISE_INTEGRATE_delta - 2>()
                            = ((step(0, j) - (2 * i + 1) * curr) / dZ1)
                            .template head<MATSLISE_INTEGRATE_delta - 2>();
                    step(0, j) = curr;
                }
                /*
                Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = stepUp(1);
                stepUp(1)(0) = 0;
                stepUp(1)(1) = 0;
                stepUp(1).template tail<MATSLISE_INTEGRATE_delta - 2>()
                        = ((stepUp(0) - (2 * j + 1) * curr) / dZ2)
                        .template head<MATSLISE_INTEGRATE_delta - 2>();
                stepUp(0) = curr;
                 */
            }
        }

        /*
         * {
            std::cout << "\n****\n" << delta << ", " << dZ1 << ", " << dZ2 << std::endl;

            Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, MATSLISE_INTEGRATE_delta> coutI;
            for (int i = 0; i < MATSLISE_ETA_delta; ++i)
                for (int j = 0; j < MATSLISE_ETA_delta; ++j)
                    coutI.row(i * MATSLISE_ETA_delta + j) = I(i, j);
            std::cout << coutI << std::endl;
            std::cout << std::endl;
        }
         */
        return I;
    } // else
    {
        Scalar I000, I101, I011, I112;

        if (equal) {
            Scalar *eta = calculateEta<Scalar>(delta * delta * dZ1, 2);

            I(0, 0)(0) = I000 = (delta * eta[0] * eta[1] + delta) / 2;
            I(0, 1)(1) = I(1, 0)(1) = I101 = I011 = (eta[0] * eta[0] - 1) / (2 * dZ1);
            I(1, 1)(2) = I112 = (delta * eta[0] * eta[1] - delta) / (2 * dZ1);

            delete[] eta;
        } else {
            Scalar denominator = 1 / (dZ1 - dZ2);
            Scalar *eta1 = calculateEta<Scalar>(Z1, 2);
            Scalar *eta2 = calculateEta<Scalar>(Z2, 2);

            I(0, 0)(0) = I000 = denominator * delta * (dZ1 * eta1[1] * eta2[0] - dZ2 * eta1[0] * eta2[1]);
            I(1, 0)(1) = I101 = denominator * (eta1[0] * eta2[0] - Z2 * eta1[1] * eta2[1] - 1);
            I(0, 1)(1) = I011 = -denominator * (eta1[0] * eta2[0] - Z1 * eta1[1] * eta2[1] - 1);
            I(1, 1)(2) = I112 = denominator * delta * (eta1[0] * eta2[1] - eta1[1] * eta2[0]);

            delete[] eta1;
            delete[] eta2;
        }

        int lastN = 100;// MATSLISE_INTEGRATE_delta - 3;
        Scalar dn = std::pow(delta, lastN);
        Scalar I00n = 0;
        Scalar I01n = 0;
        Scalar I10n = 0;
        Scalar I11n = 0;
        for (int n = lastN; n > 1; --n) {
            Scalar x1 = I101 * dn - I10n;
            Scalar x2 = I011 * dn - I01n;
            Scalar y1 = I000 * dn - I00n;
            Scalar y2 = I112 * dn - I11n;

            I00n = (x1 * dZ1 + x2 * dZ2 + dn) / n;
            I11n = (x1 + x2) / n;
            I01n = (y1 + y2 * dZ1) / n;
            I10n = (y1 + y2 * dZ2) / n;

            if (n - 1 < MATSLISE_INTEGRATE_delta)
                I(0, 0)(n - 1) = I00n;
            if (n + 1 < MATSLISE_INTEGRATE_delta)
                I(1, 1)(n + 1) = I11n;
            if (n < MATSLISE_INTEGRATE_delta) {
                I(0, 1)(n) = I01n;
                I(1, 0)(n) = I10n;
            }

            dn /= delta;
        }
    }


    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        if (i > 1) {
            for (int j = 0; j < 2; ++j)
                I(i, j).template bottomRows<MATSLISE_INTEGRATE_delta - 2>()
                        = ((I(i - 2, j) - (2 * i - 3) * I(i - 1, j)) / dZ1)
                        .template topRows<MATSLISE_INTEGRATE_delta - 2>();
        }
        for (int j = 2; j < MATSLISE_ETA_delta; ++j) {
            I(i, j).template bottomRows<MATSLISE_INTEGRATE_delta - 2>()
                    = ((I(i, j - 2) - (2 * j - 3) * I(i, j - 1)) / dZ2)
                    .template topRows<MATSLISE_INTEGRATE_delta - 2>();
        }
    }


/*
    {
        std::cout << "\n****\n" << delta << ", " << dZ1 << ", " << dZ2 << std::endl;

        Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, MATSLISE_INTEGRATE_delta> coutI;
        for (int i = 0; i < MATSLISE_ETA_delta; ++i)
            for (int j = 0; j < MATSLISE_ETA_delta; ++j)
                coutI.row(i * MATSLISE_ETA_delta + j) = I(i, j);
        std::cout << coutI << std::endl;
        std::cout << std::endl;
    }
    */

    return I;
}

template<typename Scalar, int N = MATSLISE_N>
void legendreTransform(const Scalar &delta, Eigen::Array<Scalar, N, 1> &quadratures) {
    static Eigen::Matrix<Scalar, N, N> legendrePolynomials = (
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
    ).finished().template topLeftCorner<N, N>();

    Scalar d = 1;
    for (int i = 0; i < N; ++i, d /= delta)
        quadratures(i) *= d;
    quadratures = (legendrePolynomials * quadratures.matrix()).array();
}

template<typename Scalar, bool equal = false>
Eigen::Array<Scalar, MATSLISE_N, 1> vbar_formulas(
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &y1,
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &y2,
        const Scalar &delta, const Scalar &dZ1, const Scalar &dZ2) {

    Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta>
            eta = eta_integrals<Scalar, equal>(delta, dZ1, dZ2);

    Eigen::Array<Scalar, MATSLISE_N, 1> quadratures = Eigen::Array<Scalar, MATSLISE_N, 1>::Zero();
/*
    std::cout << "\ny1: " << std::endl;
    std::cout << y1 << std::endl;
    std::cout << "\ny2: " << std::endl;
    std::cout << y2 << std::endl;
*/

    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        for (int ni = (i == 0 ? 0 : 2 * i - 1); ni < (i == 0 ? 1 : MATSLISE_HMAX_delta); ++ni) {
            for (int j = 0; j < MATSLISE_ETA_delta; ++j) {
                for (int nj = (j == 0 ? 0 : 2 * j - 1);
                     (j == 0 ? nj < 1 : nj < MATSLISE_HMAX_delta) && ni + nj < MATSLISE_INTEGRATE_delta; ++nj) {
                    Scalar s = y1(i, ni) * y2(j, nj);
                    int count = std::min(MATSLISE_N, MATSLISE_INTEGRATE_delta - ni - nj);
                    quadratures.segment(0, count) += s * eta(i, j).segment(ni + nj, count);
                }
            }
        }
    }

    //std::cout << quadratures.transpose() << std::endl;

    legendreTransform(delta, quadratures);
    return quadratures;
}

template<typename Scalar>
Eigen::Array<Scalar, MATSLISE_N, 1> vbar_formulas(
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &y,
        const Scalar &delta, const Scalar &dZ) {
    return vbar_formulas<Scalar, true>(y, y, delta, dZ, dZ);
}

#endif //MATSLISE_SECTOR_CORRECTION_POTENTIAL_H