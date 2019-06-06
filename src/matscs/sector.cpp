#include <functional>
#include "../matslise.h"
#include "matscs_formulas.h"
#include "../util/legendre.h"
#include "../util/calculateEta.h"
#include "../util/horner.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;
using namespace Eigen;

template<typename Scalar>
Matscs<Scalar>::Sector::Sector(const Matscs *s, Scalar min, Scalar max, bool backward)
        : s(s), min(min), max(max), backward(backward) {
    h = max - min;
    vs = legendre::getCoefficients(MATSCS_N, s->V, min, max);
    SelfAdjointEigenSolver<Matrix<Scalar, Dynamic, Dynamic>> es(vs[0]);
    diagonalize = es.eigenvectors();

    for (int i = 0; i < MATSCS_N; ++i)
        vs[i] = (backward && i % 2 == 1 ? -1 : 1) * diagonalize.transpose() * vs[i] * diagonalize;

    calculateTCoeffs();
}

template<typename Scalar>
void Matscs<Scalar>::Sector::calculateTCoeffs() {
    calculate_tcoeff_matrix(s->n, h, vs, t_coeff, t_coeff_h);
}

template<typename Scalar>
T<Scalar, Dynamic> Matscs<Scalar>::Sector::calculateT(Scalar E, Scalar delta, bool use_h) const {
    Matrix<Scalar, Dynamic, Dynamic> zero = Matrix<Scalar, Dynamic, Dynamic>::Zero(s->n, s->n);
    Matrix<Scalar, Dynamic, Dynamic> one = Matrix<Scalar, Dynamic, Dynamic>::Identity(s->n, s->n);

    if (fmath<Scalar>::abs(delta) <= EPS) {
        return T<Scalar, Dynamic>(s->n);
    }
    if (use_h && fmath<Scalar>::abs(delta - h) <= EPS)
        return calculateT(E);

    Matrix<Scalar, Dynamic, 1> VEd = (vs[0].diagonal() - Matrix<Scalar, Dynamic, 1>::Constant(s->n, E)) * delta;
    Matrix<Scalar, Dynamic, Dynamic> *eta = calculateEta<Scalar>(VEd * delta, s->n, MATSCS_ETA_delta);
    T<Scalar, Dynamic> t(s->n);
    t.t << zero, zero, eta[1] * VEd.asDiagonal(), zero;
    t.dt << zero, zero, -delta * eta[1] - (delta * delta / 2) * eta[2] * VEd.asDiagonal(), zero;

    for (int i = 0; i < MATSCS_ETA_delta; ++i) {
        Matrix<Scalar, Dynamic, Dynamic> hor = horner<Matrix<Scalar, Dynamic, Dynamic>>(t_coeff.row(i), delta, MATSCS_HMAX_delta);
        t.t += hor * kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), eta[i]);

        if (i + 1 < MATSCS_ETA_delta)
            t.dt += hor * kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), eta[i + 1] * (-delta * delta / 2.));
    }


    delete[] eta;

    t.t = kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize) * t.t *
          kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize.transpose());
    t.dt = kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize) * t.dt *
           kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize.transpose());
    return t;
}

template<typename Scalar>
T<Scalar, Dynamic> Matscs<Scalar>::Sector::calculateT(Scalar E, bool use_h) const {
    if (!use_h)
        return calculateT(E, h, false);
    int N = s->n;
    Matrix<Scalar, Dynamic, Dynamic> zero = Matrix<Scalar, Dynamic, Dynamic>::Zero(N, N);
    Matrix<Scalar, Dynamic, Dynamic> one = Matrix<Scalar, Dynamic, Dynamic>::Identity(N, N);

    Matrix<Scalar, Dynamic, 1> VEd = (vs[0].diagonal() - Matrix<Scalar, Dynamic, 1>::Constant(N, E)) * h;
    Matrix<Scalar, Dynamic, Dynamic> *eta = calculateEta<Scalar>(VEd * h, N, MATSCS_ETA_h);
    T<Scalar, Dynamic> t(N);
    t.t << zero, zero, eta[1] * VEd.asDiagonal(), zero;
    t.dt << zero, zero, -h * eta[1] - (h * h / 2) * eta[2] * VEd.asDiagonal(), zero;

    for (int i = 0; i < MATSCS_ETA_h; ++i) {
        t.t += t_coeff_h[i] * kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), eta[i]);

        if (i + 1 < MATSCS_ETA_h)
            t.dt += t_coeff_h[i] * kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), eta[i + 1] * (-h * h / 2));
    }
    delete[] eta;

    t.t = kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize) * t.t *
          kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize.transpose());
    t.dt = kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize) * t.dt *
           kroneckerProduct(Matrix<Scalar, 2, 2>::Identity(), diagonalize.transpose());
    return t;
}

template<typename Scalar, int r>
Y<Scalar, Dynamic, r>
propagate_delta(const typename Matscs<Scalar>::Sector *ms, Scalar E, const Y<Scalar, Dynamic, r> &y0, Scalar delta,
                bool use_h) {
    if (ms->backward)
        delta *= -1;
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    if (delta > ms->h)
        delta = ms->h;

    const T<Scalar, Dynamic> &t = ms->calculateT(E, delta, use_h);
    Y<Scalar, Dynamic, r> y = y0;
    if (ms->backward)
        y.reverse();
    Y<Scalar, Dynamic, r> y1 = forward ? t * y : t / y;
    if (ms->backward)
        y1.reverse();

    return y1;
}

template<typename Scalar>
template<int r>
Y<Scalar, Dynamic, r> Matscs<Scalar>::Sector::propagate(
        Scalar E, const Y<Scalar, Dynamic, r> &y0, Scalar a, Scalar b, bool use_h) const {
    Y<Scalar, Dynamic, r> y = y0;
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


template Y<double, Dynamic, -1>
Matscs<double>::Sector::propagate<-1>(double E, const Y<double, Dynamic, -1> &y0, double a, double b, bool use_h) const;

template Y<double, Dynamic, 1>
Matscs<double>::Sector::propagate<1>(double E, const Y<double, Dynamic, 1> &y0, double a, double b, bool use_h) const;


template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic>
propagatePsi_delta(const typename Matscs<Scalar>::Sector *sector, Scalar E, const Matrix<Scalar, Dynamic, Dynamic> &psi,
                   Scalar delta) {
    if (sector->backward)
        delta *= -1;
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    if (delta > sector->h)
        delta = sector->h;

    // TODO: verify
    Scalar extra = sector->backward ? -1 : 1;
    if (delta > 0) {
        T<Scalar, Dynamic> T = sector->calculateT(E, delta);
        return extra * (T.getT(1, 1) + extra * T.getT(1, 0) * psi).transpose()
                .colPivHouseholderQr()
                .solve((T.getT(0, 1) + extra * T.getT(0, 0) * psi).transpose())
                .transpose();
    } else if (delta < 0) {
        T<Scalar, Dynamic> T = sector->calculateT(E, -delta);
        return extra * (T.getT(1, 1) - extra * psi * T.getT(0, 1))
                .colPivHouseholderQr()
                .solve((extra * psi * T.getT(0, 0) - T.getT(1, 0)));
    } else {
        return psi;
    }
}

template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic>
Matscs<Scalar>::Sector::propagatePsi(Scalar E, const Matrix<Scalar, Dynamic, Dynamic> &_psi, Scalar a, Scalar b) const {
    Matrix<Scalar, Dynamic, Dynamic> psi = _psi;
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

template<typename Scalar>
Scalar Matscs<Scalar>::Sector::calculateError() const {
    Scalar E = vs[0].diagonal().minCoeff();
    Scalar error = (calculateT(E, true).t - calculateT(E, false).t).cwiseAbs().mean();
    if (fmath<Scalar>::isnan(error))
        return numeric_limits<Scalar>::infinity();
    return error;
}

template<typename Scalar>
Matscs<Scalar>::Sector::~Sector() {
    delete[]vs;
}

#include "../util/instantiate.h"