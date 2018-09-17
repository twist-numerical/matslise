#include <functional>
#include <Eigen/Dense>
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
    vs = legendre::getCoefficients(16, s->V, xmin, xmax);
    SelfAdjointEigenSolver<MatrixXd> es(vs[0]);
    D = es.eigenvectors();

    for (int i = 0; i < 16; ++i)
        vs[i] = D * vs[i] * D.transpose();

    calculateTCoeffs();
}

void Sector::calculateTCoeffs() {
    calculate_tcoeff_matrix(s->n, h, vs, t_coeff);

    for (int i = 0; i < MATSCS_ETA; ++i)
        t_coeff_h[i] = horner(t_coeff[i], h, MATSCS_HMAX);
}


T<MatrixXd> Sector::calculateT(double E, double delta) const {
    MatrixXd zero = MatrixXd::Zero(s->n, s->n);
    MatrixXd one = MatrixXd::Identity(s->n, s->n);

    if (fabs(delta) <= EPS) {
        return T<MatrixXd>({one, zero, zero, one}, {zero, zero, zero, zero});
    }
    if (fabs(delta - h) <= EPS)
        return calculateT(E);

    VectorXd VEd = (vs[0].diagonal() - VectorXd::Constant(s->n, E)) * delta;
    MatrixXd *eta = calculateEta(VEd * delta, s->n, MATSCS_ETA);
    T<MatrixXd> t({zero, zero, eta[1] * VEd.asDiagonal(), zero}, {zero, zero, zero, zero});

    for (int i = 0; i < MATSCS_ETA; ++i) {
        Matrix2D<MatrixXd> hor = horner(t_coeff[i], delta, MATSCS_HMAX);
        t.t += hor * eta[i];
    }


    delete[] eta;

    t.t = D.transpose() * t.t * D;
    return t;
}

T<MatrixXd> Sector::calculateT(double E) const {
    MatrixXd zero = MatrixXd::Zero(s->n, s->n);
    MatrixXd one = MatrixXd::Identity(s->n, s->n);

    VectorXd VEd = (vs[0].diagonal() - VectorXd::Constant(s->n, E)) * h;
    MatrixXd *eta = calculateEta(VEd * h, s->n, MATSCS_ETA);
    T<MatrixXd> t({zero, zero, eta[1] * VEd.asDiagonal(), zero}, {zero, zero, zero, zero});

    for (int i = 0; i < MATSCS_ETA; ++i)
        t.t += t_coeff_h[i] * eta[i];
    delete[] eta;

    t.t = D.transpose() * t.t * D;

    return t;
}

template<typename Type, int... Args>
Y<Matrix<Type, Args...>> Sector::propagate(double E, const Y<Matrix<Type, Args...>> &y0, double delta) const {
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;

    T<MatrixXd> t = calculateT(E, delta);
    return forward ? t * y0 : t / y0;
}

template Y<Matrix<double, -1, -1>> Sector::propagate(double E, const Y<Matrix<double, -1, -1>> &y0, double delta) const;
template Y<Matrix<double, -1, 1>> Sector::propagate(double E, const Y<Matrix<double, -1, 1>> &y0, double delta) const;

Sector::~Sector() {
    delete[]vs;
}
