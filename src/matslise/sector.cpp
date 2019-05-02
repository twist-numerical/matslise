#define _USE_MATH_DEFINES

#include <cmath>
#include <array>
#include "matslise_formulas.h"
#include "../util/legendre.h"
#include "../util/calculateEta.h"
#include "../util/horner.h"
#include "../util/theta.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;

Matslise::Sector::Sector(Matslise *s, double min, double max, bool backward)
        : s(s), min(min), max(max), backward(backward) {
    h = max - min;
    vs = legendre::getCoefficients(MATSLISE_N, s->V, min, max);
    if (backward) {
        for (int i = 1; i < MATSLISE_N; i += 2)
            vs[i] *= -1;
    }

    calculateTCoeffs();
}

void Matslise::Sector::calculateTCoeffs() {
    calculate_tcoeff_matrix(h, vs, t_coeff, t_coeff_h);
}

T<> Matslise::Sector::calculateT(double E, double delta, bool use_h) const {
    if (fabs(delta) <= EPS)
        return T<>();
    if (use_h && fabs(delta - h) <= EPS)
        return calculateT(E);

    double *eta = calculateEta((vs[0] - E) * delta * delta, MATSLISE_ETA_delta);
    T<> t(1);
    t.t << 0, 0, (vs[0] - E) * delta * eta[1], 0;
    t.dt << 0, 0, -delta * eta[1] + -(vs[0] - E) * delta * delta * delta * eta[2] / 2, 0;

    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        Matrix2d hor = horner(t_coeff[i], delta, MATSLISE_HMAX_delta);
        t.t += hor * eta[i];

        if (i + 1 < MATSLISE_ETA_delta)
            t.dt += hor * (-delta * delta * eta[i + 1] / 2);
    }

    delete[] eta;
    return t;
}

T<> Matslise::Sector::calculateT(double E, bool use_h) const {
    if (!use_h)
        return calculateT(E, h, false);
    double *eta = calculateEta((vs[0] - E) * h * h, MATSLISE_ETA_h);
    T<> t(1);
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

double Matslise::Sector::prufer(double E, double delta, const Y<> &y0, const Y<> &y1) const {
    double theta0 = theta(y0);

    double theta1 = theta(y1);
    double ff = E - vs[0];
    if (ff > 0) {
        double f = sqrt(ff);
        double C = atan_safe(y0.y[0] * f, y0.y[1]) / f;
        theta0 -= round(((C + delta) * f - theta1) / M_PI) * M_PI;
    } else {
        double s0 = y0.y[0] * y0.y[1];
        double s1 = y1.y[0] * y1.y[1];
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

Y<> propagate_delta(const Matslise::Sector *ms, double E, const Y<> &y0, double delta, double &theta, bool use_h) {
    if (ms->backward)
        delta *= -1;
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    if (delta > ms->h)
        delta = ms->h;

    const T<> &t = ms->calculateT(E, delta, use_h);
    Y<> y = y0;
    if (ms->backward)
        y.reverse();
    Y<> y1 = forward ? t * y : t / y;
    if (ms->backward)
        y1.reverse();

    if (forward != ms->backward)
        theta += ms->prufer(E, delta, y0, y1);
    else
        theta -= ms->prufer(E, delta, y1, y0);

    return y1;
}

Y<> Matslise::Sector::propagate(double E, const Y<> &y0, bool forward, bool use_h) const {
    double theta = 0;
    return propagate_delta(this, E, y0, forward ? h : -h, theta, use_h);
}

Y<> Matslise::Sector::propagate(double E, const Y<> &y0, double a, double b, double &theta, bool use_h) const {
    Y<> y = y0;
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

double Matslise::Sector::calculateError() const {
    return (calculateT(vs[0], true).t - calculateT(vs[0], false).t).cwiseAbs().sum();
}

Matslise::Sector::~Sector() {
    delete[]vs;
}
