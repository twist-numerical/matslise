//
// Created by toon on 5/19/20.
//

#ifndef MATSLISE_ETAINTEGRALS_H
#define MATSLISE_ETAINTEGRALS_H

#include "../util/calculateEta.h"

#define MATSLISE_INTEGRATE_delta (2*MATSLISE_HMAX_delta)


inline constexpr Eigen::Index ETA_index(Eigen::Index i, Eigen::Index j) {
    Eigen::Index n = i + j;
    Eigen::Index t = n * (n + 1) / 2 + i; //i + (j * MATSLISE_ETA_delta);
    if (i + j >= MATSLISE_ETA_delta) {
        Eigen::Index k = i + j - MATSLISE_ETA_delta + 1;
        t -= k * k;
    }
    return t;
}

template<typename Scalar, Eigen::Index N>
Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, Eigen::Dynamic> etaProduct(
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &u,
        const Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &v) {

    Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, Eigen::Dynamic> uv
            = Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, Eigen::Dynamic>
            ::Zero(MATSLISE_ETA_delta * MATSLISE_ETA_delta, N);

    for (int i = 0; i < MATSLISE_ETA_delta; ++i)
        for (int j = 0; j < MATSLISE_ETA_delta; ++j) {
            Eigen::Index ij = ETA_index(i, j);
            for (int ni = 0; ni < MATSLISE_HMAX_delta; ++ni)
                for (int nj = 0; nj < MATSLISE_HMAX_delta && ni + nj < N; ++nj)
                    uv(ij, ni + nj) += u(i, ni) * v(j, nj);
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

    Eigen::Array<Scalar, 2, 1> eta = calculateEta<Scalar, 2>(delta * delta * dZ1);

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

    return I;
}

template<typename Scalar, bool equal>
Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, Eigen::Dynamic>
eta_integrals(const Scalar &delta, const Scalar &dZ1, const Scalar &dZ2) {
    Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, Eigen::Dynamic> I
            = Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, Eigen::Dynamic>
            ::Zero(MATSLISE_ETA_delta * MATSLISE_ETA_delta, MATSLISE_INTEGRATE_delta);

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
                    I.row(ETA_index(i, j)) = stepUp(0);
            }
        }
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
                    I.row(ETA_index(i, j)) = stepUp(0);
                Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1> curr = stepUp(1);
                stepUp(1)(0) = 0;
                stepUp(1)(1) = 0;
                stepUp(1).template tail<MATSLISE_INTEGRATE_delta - 2>()
                        = ((stepUp(0) - (2 * j + 1) * curr) / dZ2)
                        .template head<MATSLISE_INTEGRATE_delta - 2>();
                stepUp(0) = curr;
            }
        }

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
                    I.row(ETA_index(i, j)) = stepUp(0);
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
            }
        }

        return I;
    } else {
        Scalar I000, I101, I011, I112;

        if (equal) {
            Eigen::Array<Scalar, 2, 1> eta = calculateEta<Scalar, 2>(delta * delta * dZ1);

            I(ETA_index(0, 0), 0) = I000 = (delta * eta[0] * eta[1] + delta) / 2;
            I(ETA_index(0, 1), 1) = I(ETA_index(1, 0), 1) = I101 = I011 = (eta[0] * eta[0] - 1) / (2 * dZ1);
            I(ETA_index(1, 1), 2) = I112 = (delta * eta[0] * eta[1] - delta) / (2 * dZ1);
        } else {
            Scalar denominator = 1 / (dZ1 - dZ2);
            Eigen::Array<Scalar, 2, 1> eta1 = calculateEta<Scalar, 2>(Z1);
            Eigen::Array<Scalar, 2, 1> eta2 = calculateEta<Scalar, 2>(Z2);

            I(ETA_index(0, 0), 0) = I000 = denominator * delta * (dZ1 * eta1[1] * eta2[0] - dZ2 * eta1[0] * eta2[1]);
            I(ETA_index(1, 0), 1) = I101 = denominator * (eta1[0] * eta2[0] - Z2 * eta1[1] * eta2[1] - 1);
            I(ETA_index(0, 1), 1) = I011 = -denominator * (eta1[0] * eta2[0] - Z1 * eta1[1] * eta2[1] - 1);
            I(ETA_index(1, 1), 2) = I112 = denominator * delta * (eta1[0] * eta2[1] - eta1[1] * eta2[0]);
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
                I(ETA_index(0, 0), n - 1) = I00n;
            if (n + 1 < MATSLISE_INTEGRATE_delta)
                I(ETA_index(1, 1), n + 1) = I11n;
            if (n < MATSLISE_INTEGRATE_delta) {
                I(ETA_index(0, 1), n) = I01n;
                I(ETA_index(1, 0), n) = I10n;
            }

            dn /= delta;
        }


        for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
            if (i > 1) {
                for (int j = 0; j < 2; ++j)
                    I.row(ETA_index(i, j)).template tail<MATSLISE_INTEGRATE_delta - 2>()
                            = ((I.row(ETA_index(i - 2, j)) - (2 * i - 3) * I.row(ETA_index(i - 1, j))) / dZ1)
                            .template head<MATSLISE_INTEGRATE_delta - 2>();
            }
            for (int j = 2; j < MATSLISE_ETA_delta; ++j) {
                I.row(ETA_index(i, j)).template tail<MATSLISE_INTEGRATE_delta - 2>()
                        = ((I.row(ETA_index(i, j - 2)) - (2 * j - 3) * I.row(ETA_index(i, j - 1))) / dZ2)
                        .template head<MATSLISE_INTEGRATE_delta - 2>();
            }
        }

        return I;
    }
}

#endif //MATSLISE_ETAINTEGRALS_H