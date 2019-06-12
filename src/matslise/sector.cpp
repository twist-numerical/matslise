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

template<typename Scalar>
Matslise<Scalar>::Sector::Sector(Matslise<Scalar> *s, const Scalar &min, const Scalar &max, bool backward)
        : s(s), min(min), max(max), backward(backward) {
    h = max - min;
    vs = legendre::getCoefficients(MATSLISE_N, s->V, min, max);
    if (backward) {
        for (int i = 1; i < MATSLISE_N; i += 2)
            vs[i] *= -1;
    }

    calculateTCoeffs();
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

    Scalar *eta = calculateEta((vs[0] - E) * delta * delta, MATSLISE_ETA_delta);
    T<Scalar> t(1);
    t.t << 0, 0, (vs[0] - E) * delta * eta[1], 0;
    t.dt << 0, 0, -delta * eta[1] + -(vs[0] - E) * delta * delta * delta * eta[2] / 2, 0;

    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        Matrix<Scalar, 2, 2> hor = horner<Matrix<Scalar, 2, 2>>(t_coeff.row(i), delta, MATSLISE_HMAX_delta);
        t.t += hor * eta[i];

        if (i + 1 < MATSLISE_ETA_delta)
            t.dt += hor * (-delta * delta * eta[i + 1] / 2);
    }

    delete[] eta;
    return t;
}

template<typename Scalar>
T<Scalar> Matslise<Scalar>::Sector::calculateT(const Scalar &E, bool use_h) const {
    if (!use_h)
        return calculateT(E, h, false);
    Scalar *eta = calculateEta((vs[0] - E) * h * h, MATSLISE_ETA_h);
    T<Scalar> t(1);
    t.t << 0, 0, (vs[0] - E) * h * eta[1], 0;
    t.dt << 0, 0, -h * eta[1] + -(vs[0] - E) * h * h * h * eta[2] / 2, 0;

    for (int i = 0; i < MATSLISE_ETA_h; ++i) {
        t.t += t_coeff_h[i] * eta[i];

        if (i + 1 < MATSLISE_ETA_h)
            t.dt += t_coeff_h[i] * (-h * h * eta[i + 1] / 2);
    }
    delete[] eta;
    return t;
}

template<typename Scalar>
Scalar Matslise<Scalar>::Sector::prufer(
        const Scalar &E, const Scalar &delta, const Y<Scalar> &y0, const Y<Scalar> &y1) const {
    Scalar theta0 = theta(y0);

    Scalar theta1 = theta(y1);
    Scalar ff = E - vs[0];
    if (ff > 0) {
        Scalar f = sqrt(ff);
        Scalar C = atan_safe(y0.y[0] * f, y0.y[1]) / f;
        theta0 -= round(((C + delta) * f - theta1) / M_PI) * M_PI;
    } else {
        Scalar s0 = y0.y[0] * y0.y[1];
        Scalar s1 = y1.y[0] * y1.y[1];
        if (y0.y[0] * y1.y[0] >= 0) {
            if (s0 > 0 && s1 < 0)
                theta1 += M_PI;
            else if (s0 < 0 && s1 > 0)
                theta1 -= M_PI;
        } else if (s0 * s1 > 0) {
            theta1 += M_PI;
        }
    }

    return theta1 - theta0;
}

template<typename Scalar>
Y<Scalar> propagate_delta(const typename Matslise<Scalar>::Sector *ms, Scalar E, const Y<Scalar> &y0,
                          Scalar delta, Scalar &theta, bool use_h) {
    if (ms->backward)
        delta *= -1;
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    if (delta > ms->h)
        delta = ms->h;

    const T<Scalar> &t = ms->calculateT(E, delta, use_h);
    Y<Scalar> y = y0;
    if (ms->backward)
        y.reverse();
    Y<Scalar> y1 = forward ? t * y : t / y;
    if (ms->backward)
        y1.reverse();

    if (forward != ms->backward)
        theta += ms->prufer(E, delta, y0, y1);
    else
        theta -= ms->prufer(E, delta, y1, y0);

    return y1;
}

template<typename Scalar>
Y<Scalar> Matslise<Scalar>::Sector::propagate(const Scalar &E, const Y<Scalar> &y0, bool forward, bool use_h) const {
    Scalar theta = 0;
    return propagate_delta(this, E, y0, forward ? h : -h, theta, use_h);
}

template<typename Scalar>
Y<Scalar> Matslise<Scalar>::Sector::propagate(
        const Scalar &E, const Y<Scalar> &y0, const Scalar &a, const Scalar &b, Scalar &theta, bool use_h) const {
    Y<Scalar> y = y0;
    if (!((a >= max && b >= max) || (a <= min && b <= min))) {
        if (!backward) { // forward
            if (a > min)
                y = propagate_delta(this, E, y, min - a, theta, use_h);
            if (b > min)
                y = propagate_delta(this, E, y, b - min, theta, use_h);
        } else {
            if (a < max)
                y = propagate_delta(this, E, y, max - a, theta, use_h);
            if (b < max)
                y = propagate_delta(this, E, y, b - max, theta, use_h);
        }
    }
    return y;
}

template<typename Scalar>
Scalar Matslise<Scalar>::Sector::calculateError() const {
    return (calculateT(vs[0], true).t - calculateT(vs[0], false).t).cwiseAbs().sum();
}

template<typename Scalar>
Matslise<Scalar>::Sector::~Sector() {
    delete[]vs;
}

#include "../util/instantiate.h"