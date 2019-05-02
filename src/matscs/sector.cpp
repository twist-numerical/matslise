#include <functional>
#include "matscs_formulas.h"
#include "../util/legendre.h"
#include "../util/calculateEta.h"
#include "../util/horner.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;

Matscs::Sector::Sector(const Matscs *s, double min, double max, bool backward) : s(s), min(min), max(max),
                                                                                 backward(backward) {
    h = max - min;
    vs = legendre::getCoefficients(MATSCS_N, s->V, min, max);
    SelfAdjointEigenSolver<MatrixXd> es(vs[0]);
    D = es.eigenvectors();

    for (int i = 0; i < MATSCS_N; ++i)
        vs[i] = (backward && i % 2 == 1 ? -1 : 1) * D.transpose() * vs[i] * D;

    calculateTCoeffs();
}

void Matscs::Sector::calculateTCoeffs() {
    calculate_tcoeff_matrix(s->n, h, vs, t_coeff, t_coeff_h);
}


T<Dynamic> Matscs::Sector::calculateT(double E, double delta, bool use_h) const {
    MatrixXd zero = MatrixXd::Zero(s->n, s->n);
    MatrixXd one = MatrixXd::Identity(s->n, s->n);

    if (fabs(delta) <= EPS) {
        return T<Dynamic>(s->n);
    }
    if (use_h && abs(delta - h) <= EPS)
        return calculateT(E);

    VectorXd VEd = (vs[0].diagonal() - VectorXd::Constant(s->n, E)) * delta;
    MatrixXd *eta = calculateEta(VEd * delta, s->n, MATSCS_ETA_delta);
    T<Dynamic> t(s->n);
    t.t << zero, zero, eta[1] * VEd.asDiagonal(), zero;
    t.dt << zero, zero, -delta * eta[1] - (delta * delta / 2) * eta[2] * VEd.asDiagonal(), zero;

    for (int i = 0; i < MATSCS_ETA_delta; ++i) {
        MatrixXd hor = horner(t_coeff[i], delta, MATSCS_HMAX_delta);
        t.t += hor * kroneckerProduct(Matrix2d::Identity(), eta[i]);

        if (i + 1 < MATSCS_ETA_delta)
            t.dt += hor * kroneckerProduct(Matrix2d::Identity(), eta[i + 1] * (-delta * delta / 2.));
    }


    delete[] eta;

    t.t = kroneckerProduct(Matrix2d::Identity(), D) * t.t * kroneckerProduct(Matrix2d::Identity(), D.transpose());
    t.dt = kroneckerProduct(Matrix2d::Identity(), D) * t.dt * kroneckerProduct(Matrix2d::Identity(), D.transpose());
    return t;
}

T<Dynamic> Matscs::Sector::calculateT(double E, bool use_h) const {
    if (!use_h)
        return calculateT(E, h, false);
    int N = s->n;
    MatrixXd zero = MatrixXd::Zero(N, N);
    MatrixXd one = MatrixXd::Identity(N, N);

    VectorXd VEd = (vs[0].diagonal() - VectorXd::Constant(N, E)) * h;
    MatrixXd *eta = calculateEta(VEd * h, N, MATSCS_ETA_h);
    T<Dynamic> t(N);
    t.t << zero, zero, eta[1] * VEd.asDiagonal(), zero;
    t.dt << zero, zero, -h * eta[1] - (h * h / 2) * eta[2] * VEd.asDiagonal(), zero;

    for (int i = 0; i < MATSCS_ETA_h; ++i) {
        t.t += t_coeff_h[i] * kroneckerProduct(Matrix2d::Identity(), eta[i]);

        if (i + 1 < MATSCS_ETA_h)
            t.dt += t_coeff_h[i] * kroneckerProduct(Matrix2d::Identity(), eta[i + 1] * (-h * h / 2));
    }
    delete[] eta;

    t.t = kroneckerProduct(Matrix2d::Identity(), D) * t.t * kroneckerProduct(Matrix2d::Identity(), D.transpose());
    t.dt = kroneckerProduct(Matrix2d::Identity(), D) * t.dt * kroneckerProduct(Matrix2d::Identity(), D.transpose());
    return t;
}

template<int r>
Y<Dynamic, r> propagate_delta(const Matscs::Sector *ms, double E, const Y<Dynamic, r> &y0, double delta, bool use_h) {
    if (ms->backward)
        delta *= -1;
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    if (delta > ms->h)
        delta = ms->h;

    const T<Dynamic> &t = ms->calculateT(E, delta, use_h);
    Y<Dynamic, r> y = y0;
    if (ms->backward)
        y.reverse();
    Y<Dynamic, r> y1 = forward ? t * y : t / y;
    if (ms->backward)
        y1.reverse();

    return y1;
}

template<int r>
Y<Dynamic, r> Matscs::Sector::propagate(double E, const Y<Dynamic, r> &y0, double a, double b, bool use_h) const {
    Y<Dynamic, r> y = y0;
    if (!((a >= max && b >= max) || (a <= min && b <= min))) {
        if (!backward) { // forward
            if (a > min)
                y = propagate_delta(this, E, y, min - a, use_h);
            if (b > min)
                y = propagate_delta(this, E, y, b - min, use_h);
        } else {
            if (a < max)
                y = propagate_delta(this, E, y, max - a, use_h);
            if (b < max)
                y = propagate_delta(this, E, y, b - max, use_h);
        }
    }
    return y;
}

template Y<Dynamic, -1>
Matscs::Sector::propagate<-1>(double E, const Y<Dynamic, -1> &y0, double a, double b, bool use_h) const;

template Y<Dynamic, 1>
Matscs::Sector::propagate<1>(double E, const Y<Dynamic, 1> &y0, double a, double b, bool use_h) const;


MatrixXd propagatePsi_delta(const Matscs::Sector *sector, double E, const MatrixXd &psi, double delta) {
    if (sector->backward)
        delta *= -1;
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    if (delta > sector->h)
        delta = sector->h;

    // TODO: verify
    double extra = sector->backward ? -1 : 1;
    if (delta > 0) {
        T<Dynamic> T = sector->calculateT(E, delta);
        return extra * (T.getT(1, 1) + extra * T.getT(1, 0) * psi).transpose()
                .colPivHouseholderQr()
                .solve((T.getT(0, 1) + extra * T.getT(0, 0) * psi).transpose())
                .transpose();
    } else if (delta < 0) {
        T<Dynamic> T = sector->calculateT(E, -delta);
        return extra * (T.getT(1, 1) - extra * psi * T.getT(0, 1))
                .colPivHouseholderQr()
                .solve((extra * psi * T.getT(0, 0) - T.getT(1, 0)));
    } else {
        return psi;
    }
}

MatrixXd Matscs::Sector::propagatePsi(double E, const MatrixXd &_psi, double a, double b) const {
    MatrixXd psi = _psi;
    if (!((a >= max && b >= max) || (a <= min && b <= min))) {
        if (!backward) { // forward
            if (a > min)
                psi = propagatePsi_delta(this, E, psi, min - a);
            if (b > min)
                psi = propagatePsi_delta(this, E, psi, b - min);
        } else {
            if (a < max)
                psi = propagatePsi_delta(this, E, psi, max - a);
            if (b < max)
                psi = propagatePsi_delta(this, E, psi, b - max);
        }
    }
    return psi;
}

double Matscs::Sector::calculateError() const {
    double E = vs[0].diagonal().minCoeff();
    double error = (calculateT(E, true).t - calculateT(E, false).t).cwiseAbs().mean();
    if (isnan(error))
        return numeric_limits<double>::infinity();
    return error;
}

Matscs::Sector::~Sector() {
    delete[]vs;
}
