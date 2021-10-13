#include <cmath>
#include <array>
#include "matslise_formulas.h"
#include "../util/legendre.h"
#include "../util/calculateEta.h"
#include "../util/horner.h"
#include "../util/theta.h"
#include "../util/constants.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;
using namespace Eigen;


template<typename Scalar>
int sign(const Scalar &s) {
    return s > 0 ? 1 : s < 0 ? -1 : 0;
}

template<typename Scalar>
Matslise<Scalar>::Sector::Sector(const Matslise<Scalar> *s, const Scalar &min, const Scalar &max, Direction direction)
        : min(min), max(max), direction(direction) {
    h = max - min;
    vs = legendre::getCoefficients<MATSLISE_N>(s->potential, min, max);
    if (direction == backward) {
        for (int i = 1; i < MATSLISE_N; i += 2)
            vs[i] *= -1;
    }

    if (direction != none)
        calculateTCoeffs();
}

template<typename Scalar>
void Matslise<Scalar>::Sector::setDirection(Direction newDirection) {
    if (direction != newDirection) {
        if (direction == backward || newDirection == backward)
            for (int i = 1; i < MATSLISE_N; i += 2)
                vs[i] *= -1;

        direction = newDirection;
        if (direction != none)
            calculateTCoeffs();
    }
}

template<typename Scalar>
void Matslise<Scalar>::Sector::calculateTCoeffs() {
    calculate_tcoeff_matrix(h, vs, t_coeff, t_coeff_h);
}

template<typename Scalar>
T<Scalar> Matslise<Scalar>::Sector::calculateT(const Scalar &E, const Scalar &delta, bool use_h) const {
    if (abs(delta) <= EPS)
        return T<Scalar>();
    if (use_h && abs(delta - h) <= EPS)
        return calculateT(E);

    const Scalar dd = delta * delta;
    const Scalar vsE = vs[0] - E;
    Array<Scalar, MATSLISE_ETA_delta, 1> eta = calculateEta<Scalar, MATSLISE_ETA_delta>(vsE * dd);
    T<Scalar> t(1);
    t.t << 0, 0, vsE * delta * eta[1], 0;
    t.dt << 0, 0, -delta * eta[1] - vsE * dd * delta * eta[2] / 2, 0;

    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        Matrix<Scalar, 2, 2> hor = horner<Matrix<Scalar, 2, 2>>(t_coeff.row(i), delta, MATSLISE_HMAX_delta);
        t.t += hor * eta[i];

        if (i + 1 < MATSLISE_ETA_delta)
            t.dt += hor * (-dd * eta[i + 1] / 2);
    }

    return t;
}

template<typename Scalar>
T<Scalar> Matslise<Scalar>::Sector::calculateT(const Scalar &E, bool use_h) const {
    if (!use_h)
        return calculateT(E, h, false);
    Array<Scalar, MATSLISE_ETA_h, 1> eta = calculateEta<Scalar, MATSLISE_ETA_h>((vs[0] - E) * h * h);
    T<Scalar> t(1);
    t.t << 0, 0, (vs[0] - E) * h * eta[1], 0;
    t.dt << 0, 0, -h * eta[1] + -(vs[0] - E) * h * h * h * eta[2] / 2, 0;

    for (int i = 0; i < MATSLISE_ETA_h; ++i) {
        t.t += t_coeff_h[i] * eta[i];

        if (i + 1 < MATSLISE_ETA_h)
            t.dt += t_coeff_h[i] * (-h * h * eta[i + 1] / 2);
    }

    return t;
}

template<typename Scalar>
Scalar rescale(const Scalar &theta, const Scalar &sigma) {
    Scalar sinTheta = sin(theta);
    Scalar cosTheta = cos(theta);
    return theta + atan2((sigma - 1) * sinTheta * cosTheta, 1 + (sigma - 1) * sinTheta * sinTheta);
}

template<typename Scalar>
Scalar Matslise<Scalar>::Sector::theta0(const Scalar &E, const Y<Scalar> &y0) const {
    Scalar scaling = E - v0Match > 1 ? sqrt(E - v0Match) : 1;
    if (E > vs[0]) {
        Scalar omega = sqrt(E - vs[0]);
        Scalar theta0 = atan_safe(y0.data[0] * omega, y0.data[1]);
        return rescale(theta0, scaling / omega);
    } else {
        return atan_safe(scaling * y0.data[0], y0.data[1]);
    }
}

template<typename Scalar>
Scalar Matslise<Scalar>::Sector::prufer(
        const Scalar &E, const Scalar &delta, const Y<Scalar> &y0, const Y<Scalar> &y1) const {
    Scalar scaling = E - v0Match > 1 ? sqrt(E - v0Match) : 1;
    if (E > vs[0]) {
        // page 56 (PhD Ledoux)
        Scalar omega = sqrt(E - vs[0]);
        Scalar theta0 = atan_safe(y0.data[0] * omega, y0.data[1]);
        Scalar phi_star = atan_safe(y1.data[0] * omega, y1.data[1]);
        Scalar phi_bar = omega * delta + theta0;
        phi_bar -= floor(phi_bar / constants<Scalar>::PI) * constants<Scalar>::PI;
        Scalar theta1 = phi_star - phi_bar;
        if (theta1 < -constants<Scalar>::PI / 2)
            theta1 += constants<Scalar>::PI;
        else if (theta1 > constants<Scalar>::PI / 2)
            theta1 -= constants<Scalar>::PI;
        theta1 += theta0 + omega * delta;

        if (theta1 < theta0) {
            // theta has to be increasing
            theta1 += constants<Scalar>::PI;
        }

        return rescale(theta1, scaling / omega) - rescale(theta0, scaling / omega);
    } else {
        Scalar theta0 = atan_safe(scaling * y0.data[0], y0.data[1]);
        Scalar theta1 = atan_safe(scaling * y1.data[0], y1.data[1]);
        if (y0.data[0] * y1.data[0] >= 0) {
            int signTheta0 = theta0 == 0 ? 1 : sign(theta0);
            int signTheta1 = theta1 == 0 ? -1 : sign(theta1);
            if (signTheta0 * signTheta1 < 0)
                theta1 += signTheta0 * constants<Scalar>::PI;
        } else if (theta0 * theta1 > 0) {
            theta1 += constants<Scalar>::PI;
        }

        return theta1 - theta0;
    }
}

template<typename Scalar, bool withPrufer, int cols>
void propagateDelta(const typename Matslise<Scalar>::Sector *sector,
                    const Scalar &E, Y<Scalar, 1, cols> &y, Scalar delta, bool use_h, Scalar &theta) {
    if (sector->direction == backward)
        delta *= -1;
    bool direction = delta >= 0;
    if (!direction)
        delta = -delta;
    if (delta > sector->h)
        delta = sector->h;

    Y<Scalar, 1, cols> y0 = y;
    (void) y0; // Unused when not 'withPrufer'
    const T<Scalar> &t = sector->calculateT(E, delta, use_h);
    if (sector->direction == backward)
        y.reverse();
    y = direction ? t * y : t / y;
    if (sector->direction == backward)
        y.reverse();

    if constexpr(withPrufer) {
        static_assert(cols == 1);
        theta += direction != (sector->direction == backward) ?
                 sector->prufer(E, delta, y0, y) : -sector->prufer(E, delta, y, y0);
    }
}

template<typename Scalar>
template<bool withPrufer, int cols>
typename conditional<withPrufer, std::pair<Y<Scalar, 1, cols>, Scalar>, Y<Scalar, 1, cols>>::type
Matslise<Scalar>::Sector::propagate(
        const Scalar &E, const Y<Scalar, 1, cols> &y0, const Scalar &a, const Scalar &b, bool use_h) const {
    Y<Scalar, 1, cols> y = y0;
    Scalar argdet = 0;
    if (!((a >= max && b >= max) || (a <= min && b <= min))) {
        if (direction == forward) {
            if (a > min)
                propagateDelta<Scalar, withPrufer, cols>(this, E, y, min - a, use_h, argdet);
            if (b > min)
                propagateDelta<Scalar, withPrufer, cols>(this, E, y, b - min, use_h, argdet);
        } else {
            if (a < max)
                propagateDelta<Scalar, withPrufer, cols>(this, E, y, max - a, use_h, argdet);
            if (b < max)
                propagateDelta<Scalar, withPrufer, cols>(this, E, y, b - max, use_h, argdet);
        }
    }
    if constexpr(withPrufer)
        return {y, argdet};
    else
        return y;
}

template<typename Scalar>
Scalar Matslise<Scalar>::Sector::error() const {
    Scalar e_loc0 = (
            (calculateT(vs[0], true).t - calculateT(vs[0], false).t).array() *
            (Array<Scalar, 2, 2>() << 1, 1 / h, h, 1).finished()
    ).cwiseAbs().maxCoeff();
    if constexpr (MATSLISE_N < 16)
        return e_loc0;
    Scalar h12 = pow(h, 12);
    Scalar e_locu = 0.0024 * abs(vs[15]) * h12 * h * h * h + 2.5e-5 * vs[7] * vs[7] * h12 * h * h +
                    (1.8e-5 * abs(vs[3] * vs[10]) + 2e-5 * abs(vs[2] * vs[11]) + 6e-6 * abs(vs[5] * vs[8]) +
                     4e-6 * abs(vs[4] * vs[9])) * h12 * h;
    Scalar e_locup = (0.0002 * vs[7] * vs[7] + 0.0003 * abs(vs[5] * vs[9]) + 0.0003 * abs(vs[6] * vs[8]) +
                      0.0006 * abs(vs[1] * vs[13]) + 0.0005 * abs(vs[2] * vs[12]) + 0.0004 * abs(vs[3] * vs[11])) *
                     h12 * h * h;
    Scalar e_locv = 0.0002 * abs(vs[14]) * h12 * h * h + 1e-5 * vs[6] * vs[6] * h12;


    return std::max(e_loc0, std::max(e_locu, std::max(e_locup, e_locv)));
}

#define INSTANTIATE_PROPAGATE(Scalar, withPrufer, cols) \
template typename conditional<withPrufer, pair<Y<Scalar, 1, cols>, Scalar>, Y<Scalar, 1, cols>>::type \
Matslise<Scalar>::Sector::propagate<withPrufer, cols>( \
    const Scalar&, const Y<Scalar, 1, cols>&, const Scalar&, const Scalar&, bool) const;

#define INSTANTIATE_MORE(Scalar) \
INSTANTIATE_PROPAGATE(Scalar, true, 1) \
INSTANTIATE_PROPAGATE(Scalar, false, 1) \
INSTANTIATE_PROPAGATE(Scalar, false, 2)

#include "instantiate.h"