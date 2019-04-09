#include <functional>
#include "matscs_formulas.h"
#include "../util/legendre.h"
#include "../util/calculateEta.h"
#include "../util/horner.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;
using namespace matslise::matscs_util;

Sector::Sector(const Matscs *s, double xmin, double xmax) : s(s), xmin(xmin), xmax(xmax) {
    h = xmax - xmin;
    vs = legendre::getCoefficients(MATSCS_N, s->V, xmin, xmax);
    SelfAdjointEigenSolver<MatrixXd> es(vs[0]);
    D = es.eigenvectors();

    for (int i = 0; i < MATSCS_N; ++i)
        vs[i] = D.transpose() * vs[i] * D;

    calculateTCoeffs();
}

void Sector::calculateTCoeffs() {
    calculate_tcoeff_matrix(s->n, h, vs, t_coeff, t_coeff_h);
}


T<Dynamic> Sector::calculateT(double E, double delta) const {
    MatrixXd zero = MatrixXd::Zero(s->n, s->n);
    MatrixXd one = MatrixXd::Identity(s->n, s->n);

    if (fabs(delta) <= EPS) {
        return T<Dynamic>(s->n);
    }
    if (fabs(delta - h) <= EPS)
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
    return t;
}

T<Dynamic> Sector::calculateT(double E) const {
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
    return t;
}

template<int r>
Y<Dynamic, r> Sector::propagate(double E, const Y<Dynamic, r> &y0, double delta) const {
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;

    T<Dynamic> t = calculateT(E, delta);
    return forward ? t * y0 : t / y0;
}

template Y<Dynamic, -1>
Sector::propagate<-1>(double E, const Y<Dynamic, -1> &y0, double delta) const;

template Y<Dynamic, 1>
Sector::propagate<1>(double E, const Y<Dynamic, 1> &y0, double delta) const;

MatrixXd Sector::propagatePsi(double E, const MatrixXd &psi, double delta) const {
    if (delta > 0) {
        T<Dynamic> T = calculateT(E, delta);
        return (T.getT(1, 1) + T.getT(1, 0) * psi).transpose()
                .colPivHouseholderQr()
                .solve((T.getT(0, 1) + T.getT(0, 0) * psi).transpose())
                .transpose();
    } else if (delta < 0) {
        T<Dynamic> T = calculateT(E, -delta);
        return (T.getT(0, 0) - T.getT(1, 0) * psi).transpose()
                .colPivHouseholderQr()
                .solve((-T.getT(0, 1) + T.getT(1, 1) * psi).transpose())
                .transpose();
    } else {
        return psi;
    }
}

Sector::~Sector() {
    delete[]vs;
}
