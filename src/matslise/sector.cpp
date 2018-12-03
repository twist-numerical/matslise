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
using namespace matslise::matslise_util;

Sector::Sector(Matslise *s, double xmin, double xmax) : s(s), xmin(xmin), xmax(xmax) {
    h = xmax - xmin;
    vs = legendre::getCoefficients(MATSLISE_N, s->V, xmin, xmax);

    calculateTCoeffs();
}

void Sector::calculateTCoeffs() {
    calculate_tcoeff_matrix(h, vs, t_coeff, t_coeff_h);
}

T<double> Sector::calculateT(double E, double delta) const {
    if (fabs(delta) <= EPS)
        return T<double>();
    if (fabs(delta - h) <= EPS)
        return calculateT(E);

    double *eta = calculateEta((vs[0] - E) * delta * delta, MATSLISE_ETA_delta);
    T<double> t({0, 0, (vs[0] - E) * delta * eta[1], 0},
                {0, 0, -delta * eta[1] + -(vs[0] - E) * delta * delta * delta * eta[2] / 2, 0});

    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        Matrix2D<> hor = horner(t_coeff[i], delta, MATSLISE_HMAX_delta);
        t.t += hor * eta[i];

        if (i + 1 < MATSLISE_ETA_delta)
            t.dt += hor * (-delta * delta * eta[i + 1] / 2);
    }

    delete[] eta;
    return t;
}

T<double> Sector::calculateT(double E) const {
    double *eta = calculateEta((vs[0] - E) * h * h, MATSLISE_ETA_h);
    T<double> t({0, 0, (vs[0] - E) * h * eta[1], 0},
                {0, 0, -h * eta[1] + -(vs[0] - E) * h * h * h * eta[2] / 2, 0});

    for (int i = 0; i < MATSLISE_ETA_h; ++i) {
        t.t += t_coeff_h[i] * eta[i];

        if (i + 1 < MATSLISE_ETA_h)
            t.dt += t_coeff_h[i] * (-h * h * eta[i + 1] / 2);
    }
    delete[] eta;
    return t;
}

double Sector::prufer(double E, double delta, const Y<double> &y0, const Y<double> &y1) const {
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

Y<double> Sector::propagate(double E, const Y<double> &y0, bool forward) const {
    const T<double> &t = calculateT(E);
    return forward ? t * y0 : t / y0;
}

Y<double> Sector::propagate(double E, const Y<double> &y0, double delta) const {
    if (delta >= 0)
        return calculateT(E, delta) * y0;
    else
        return calculateT(E, -delta) / y0;
}

Y<double> Sector::propagate(double E, const Y<double> &y0, double delta, double &theta) const {
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;

    const T<double> &t = calculateT(E, delta);
    Y<double> y1 = forward ? t * y0 : t / y0;

    if (forward)
        theta += prufer(E, delta, y0, y1);
    else
        theta -= prufer(E, delta, y1, y0);

    /*
    if(forward)
        cout << "(" << xmin << ")   " << xmin << " -> " << (xmin+delta) << endl;
    else
        cout << "(" << xmin << ")   " << (xmin+delta) << " -> " << xmin << endl;

    cout << "| " << y0.y[0] << ", " << y0.dy[0] << " |    | " << y1.y[0] << ", " << y1.dy[0] << " |" << endl;
    cout << "| " << y0.y[1] << ", " << y0.dy[1] << " | -> | " << y1.y[1] << ", " << y1.dy[1] << " |" << endl;
    */

    return y1;
}

Sector::~Sector() {
    delete[]vs;
}
