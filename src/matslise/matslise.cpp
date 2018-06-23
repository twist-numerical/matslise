//
// Created by toon on 5/16/18.
//

#include <cmath>
#include <array>
#include <vector>
#include <queue>
#include <Eigen/Dense>
#include "matslise_formulas.h"
#include "../matslise.h"
#include "../util/legendre.h"
#include "../util/calculateEta.h"

#define EPS (1.e-12)

using namespace matslise;
using namespace std;
using namespace Eigen;

Matslise::Matslise(function<double(double)> V, double xmin, double xmax, int sectorCount)
        : V(V), xmin(xmin), xmax(xmax), sectorCount(sectorCount) {
    sectors = new Sector *[sectorCount];
    double h = (xmax - xmin) / sectorCount;
    double mid = (xmax + xmin) / 2;
    for (int i = 0; i < sectorCount; ++i) {
        double a = xmin + i * h;
        double b = a + h;
        if (b - 1.e-5 > mid) {
            match = b;
            mid = xmax + 1;
        }
        sectors[i] = new Sector(this, a, b);
    }
}


tuple<Y, double> Matslise::propagate(double E, const Y &_y, double a, double b) const {
    Y y = _y;
    double theta = y.theta();
    if (a < b) {
        for (int i = 0; i < sectorCount; ++i) {
            Sector *sector = sectors[i];
            if (sector->xmax > a) {
                if (sector->xmin < a) // first
                    y = sector->propagate(E, y, sector->xmin - a, theta);

                if (sector->xmax > b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, theta);
                    break;
                }

                y = sector->propagate(E, y, sector->xmax - sector->xmin, theta);
            }
        }
    } else {
        for (int i = sectorCount - 1; i >= 0; --i) {
            Sector *sector = sectors[i];
            if (sector->xmin < a) {
                if (sector->xmax > a) // first
                    y = sector->propagate(E, y, sector->xmin - a, theta);
                else
                    y = sector->propagate(E, y, sector->xmin - sector->xmax, theta);

                if (sector->xmin < b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, theta);
                    break;
                }

            }
        }
    }
    return make_tuple(y, theta);
}

tuple<double, double, double>
Matslise::calculateError(double E, const Y &left, const Y &right) const {
    Y l, r;
    double thetaL, thetaR;
    tie(l, thetaL) = propagate(E, left, xmin, match);
    tie(r, thetaR) = propagate(E, right, xmax, match);
    return make_tuple(l.y[1] * r.y[0] - r.y[1] * l.y[0],
                      l.dy[1] * r.y[0] + l.y[1] * r.dy[0] - (r.dy[1] * l.y[0] + r.y[1] * l.dy[0]),
                      thetaL - thetaR);
}

tuple<unsigned int, double> newtonIteration(const Matslise *ms, double E, const Y &left, const Y &right, double tol) {
    double adjust, error, derror, theta;
    int i = 0;
    do {
        tie(error, derror, theta) = ms->calculateError(E, left, right);
        adjust = error / derror;
        E = E - adjust;
        if (++i > 50) {
            throw runtime_error("Newton-iteration did not converge");
        }
    } while (fabs(adjust) > tol);

    int index = (int) (round(theta / M_PI) - 1);
    if (index < 0)
        index = 0;
    return make_tuple(index, E);
}

vector<tuple<unsigned int, double>> *
Matslise::computeEigenvaluesByIndex(unsigned int Imin, unsigned int Imax, const Y &left, const Y &right) const {
    double Emin = -1;
    double Emax = 1;
    while (true) {
        unsigned int i = (unsigned int) floor(get<2>(calculateError(Emax, left, right)) / M_PI);
        if (i >= Imax)
            break;
        else {
            if (i < Imin)
                Emin = Emax;
            Emax *= 2;
        }
    }
    if (Emin == -1) {
        while (true) {
            unsigned int i = (unsigned int) floor(get<2>(calculateError(Emin, left, right)) / M_PI);
            if (i <= Imin)
                break;
            else {
                if (i > Imax)
                    Emax = Emin;
                Emin *= 2;
            }
        }
    }
    return computeEigenvalues(Emin, Emax, Imin, Imax, left, right);
};

vector<tuple<unsigned int, double>> *
Matslise::computeEigenvalues(double Emin, double Emax, const Y &left, const Y &right) const {
    return computeEigenvalues(Emin, Emax, 0, UINT_MAX, left, right);
};

vector<tuple<unsigned int, double>> *
Matslise::computeEigenvalues(double Emin, double Emax, unsigned int Imin, unsigned int Imax, const Y &left,
                             const Y &right) const {
    vector<tuple<unsigned int, double>> *eigenvalues = new vector<tuple<unsigned int, double>>();
    queue<tuple<double, double, double, double, unsigned int>> toCheck;

    toCheck.push(make_tuple(Emin, get<2>(calculateError(Emin, left, right)) / M_PI,
                            Emax, get<2>(calculateError(Emax, left, right)) / M_PI,
                            0));

    double a, ta, b, tb, c, tc;
    unsigned int ia, ib, depth;
    while (!toCheck.empty()) {
        tie(a, ta, b, tb, depth) = toCheck.front();
        toCheck.pop();
        ia = (unsigned int) floor(ta);
        ib = (unsigned int) floor(tb);
        if (ta >= tb || ia == ib || ib <= Imin || Imax <= ia)
            continue;
        if (ia + 1 == ib)
            ++depth;

        c = (a + b) / 2;
        if (tb - ta < 0.1 || depth > 10)
            eigenvalues->push_back(newtonIteration(this, c, left, right, 1e-12));
        else {
            tc = get<2>(calculateError(c, left, right)) / M_PI;
            toCheck.push(make_tuple(a, ta, c, tc, depth));
            toCheck.push(make_tuple(c, tc, b, tb, depth));
        }
    }

    sort(eigenvalues->begin(), eigenvalues->end());

    return eigenvalues;
}

Matslise::~Matslise() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
}

Array<Y, Dynamic, 1> Matslise::computeEigenfunction(double E, const matslise::Y &left, const matslise::Y &right,
                                                    const ArrayXd &x) const {
    long n = x.size();
    for (int i = 1; i < n; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");

    Array<Y, Dynamic, 1> ys(n);

    long forward = 0;
    Y y = left;
    int iLeft = 0;
    { // left
        while (forward < n && x[forward] < xmin - EPS)
            ++forward;

        Sector *sector = sectors[0];
        for (; forward < n; ++forward) {
            while (sector->xmax < x[forward]) {
                y = sector->calculateT(E) * y;
                ++iLeft;
                sector = sectors[iLeft];
                if (iLeft >= sectorCount || sector->xmin > match)
                    goto allLeftSectorsDone;
            }

            ys[forward] = sector->calculateT(E, x[forward] - sector->xmin) * y;
        }
        allLeftSectorsDone:;
    }

    { // right
        Y yLeft = y;

        long reverse = n;
        while (reverse-- > forward && x[reverse] > xmax + EPS);
        long lastValid = reverse;

        Sector *sector = sectors[sectorCount - 1];
        y = sector->calculateT(E) / right;
        for (int i = sectorCount - 1; reverse >= 0; --reverse) {
            while (sector->xmin > x[reverse]) {
                --i;
                sector = sectors[i];
                if (i < iLeft || sector->xmax < match)
                    goto allRightSectorsDone;
                y = sector->calculateT(E) / y;
            }

            ys[reverse] = sector->calculateT(E, x[reverse] - sector->xmin) * y;
        }
        allRightSectorsDone:;

        Y yRight = y;
        double scale = yLeft.y[0] / yRight.y[0];
        for (long i = reverse + 1; i <= lastValid; ++i) {
            ys[i].y *= scale;
            ys[i].dy *= scale;
        }
    }

    return ys;
}

Sector::Sector(Matslise *s, double xmin, double xmax) : s(s), xmin(xmin), xmax(xmax) {
    h = xmax - xmin;
    vs = legendre::getCoefficients(17, s->V, xmin, xmax);

    calculateTCoeffs();
}

void Sector::calculateTCoeffs() {
    double v1 = vs[1],
            v2 = vs[2],
            v3 = vs[3],
            v4 = vs[4],
            v5 = vs[5],
            v6 = vs[6],
            v7 = vs[7],
            v8 = vs[8],
            v9 = vs[9],
            v10 = vs[10],
            v11 = vs[11],
            v12 = vs[12],
            v13 = vs[13],
            v14 = vs[14],
            v15 = vs[15],
            v16 = vs[16];

    // @formatter:off
    u = MATSLISE_U;
    up = MATSLISE_UP;
    v = MATSLISE_V;
    vp = MATSLISE_VP;
    // @formatter:on

    for (int i = 0; i < MATSLISE_ETA; ++i) {
        hu[i] = hup[i] = hv[i] = hvp[i] = 0;
        double H = 1;
        for (int j = 0; j < MATSLISE_HMAX; ++j, H *= h) {
            hu[i] += H * u[i][j];
            hup[i] += H * up[i][j];
            hv[i] += H * v[i][j];
            hvp[i] += H * vp[i][j];
        }
    }
}

T Sector::calculateT(double E, double delta) const {
    if (fabs(delta) <= EPS)
        return T();
    if (fabs(delta - h) <= EPS)
        return calculateT(E);

    double *eta = calculateEta((vs[0] - E) * delta * delta, MATSLISE_ETA);
    T t((Matrix2d() << 0, 0, (vs[0] - E) * delta * eta[1], 0).finished(),
        (Matrix2d() << 0, 0, -delta * eta[1] + -(vs[0] - E) * delta * delta * delta * eta[2] / 2, 0).finished());

    for (int i = 0; i < MATSLISE_ETA; ++i) {
        double D = 1;
        for (int j = 0; j < MATSLISE_HMAX; ++j, D *= delta) {
            t.t(0, 0) += D * eta[i] * u[i][j];
            t.t(0, 1) += D * eta[i] * v[i][j];
            t.t(1, 0) += D * eta[i] * up[i][j];
            t.t(1, 1) += D * eta[i] * vp[i][j];

            if (i + 1 < MATSLISE_ETA) {
                double dEta = -delta * delta * eta[i + 1] / 2;
                t.dt(0, 0) += D * dEta * u[i][j];
                t.dt(0, 1) += D * dEta * v[i][j];
                t.dt(1, 0) += D * dEta * up[i][j];
                t.dt(1, 1) += D * dEta * vp[i][j];
            }
        }
    }

    delete[] eta;
    return t;
}

T Sector::calculateT(double E) const {
    double *eta = calculateEta((vs[0] - E) * h * h, MATSLISE_ETA);
    T t((Matrix2d() << 0, 0, (vs[0] - E) * h * eta[1], 0).finished(),
        (Matrix2d() << 0, 0, -h * eta[1] + -(vs[0] - E) * h * h * h * eta[2] / 2, 0).finished());

    for (int i = 0; i < MATSLISE_ETA; ++i) {
        t.t(0, 0) += eta[i] * hu[i];
        t.t(0, 1) += eta[i] * hv[i];
        t.t(1, 0) += eta[i] * hup[i];
        t.t(1, 1) += eta[i] * hvp[i];

        if (i + 1 < MATSLISE_ETA) {
            double dEta = -h * h * eta[i + 1] / 2;
            t.dt(0, 0) += dEta * hu[i];
            t.dt(0, 1) += dEta * hv[i];
            t.dt(1, 0) += dEta * hup[i];
            t.dt(1, 1) += dEta * hvp[i];
        }
    }
    delete[] eta;
    return t;
}

double Sector::prufer(double E, double delta, const Y &y0, const Y &y1) const {
    double theta0 = y0.theta();

    double theta1 = y1.theta();
    double ff = E - vs[0];
    if (ff > 0) {
        double f = sqrt(ff);
        double C = atan(y0.y[0] / y0.y[1] * f) / f;
        theta1 += round(((C + delta) * f - theta1) / M_PI) * M_PI;
    } else {
        if (y0.y[0] * y1.y[0] >= 0) {
            if (theta0 > 0 && theta1 < 0)
                theta1 += M_PI;
            else if (theta0 < 0 && theta1 > 0)
                theta1 -= M_PI;
        } else if (theta0 * theta1 > 0) {
            theta1 += M_PI;
        }
    }

    return theta1 - theta0;
}

Y Sector::propagate(double E, const Y &y0, double delta, double &theta) const {
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    const T &t = calculateT(E, delta);
    const Y y1 = forward ? t * y0 : t / y0;

    if (forward)
        theta += prufer(E, delta, y0, y1);
    else
        theta -= prufer(E, delta, y1, y0);

    return y1;
}

Sector::~Sector() {
    delete[]vs;
}
