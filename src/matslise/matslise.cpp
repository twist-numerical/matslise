//
// Created by toon on 5/16/18.
//

#include <array>
#include <vector>
#include <queue>
#include "../matslise.h"
#include "matslise_formulas.h"
#include "../util/theta.h"

#define EPS (1.e-12)

using namespace matslise;
using namespace matslise::matslise_util;
using namespace std;
using namespace Eigen;

Matslise::Matslise(function<double(double)> V, double xmin, double xmax, int sectorCount)
        : V(V), xmin(xmin), xmax(xmax), sectorCount(sectorCount) {
    sectors = new Sector *[sectorCount];
    double h = (xmax - xmin) / sectorCount;
    double mid = .51 * xmax + .49 * xmin;
    for (int i = 0; i < sectorCount; ++i) {
        double a = xmin + i * h;
        double b = xmax - (sectorCount - i - 1) * h;
        if (b > mid) {
            match = b;
            mid = xmax + 1;
        }
        sectors[i] = new Sector(this, a, b);
    }
}


pair<Y<>, double> Matslise::propagate(double E, const Y<> &_y, double a, double b) const {
    Y<> y = _y;
    double theta = matslise::theta(y);
    if (a <= b) {
        for (int i = 0; i < sectorCount; ++i) {
            Sector *sector = sectors[i];
            if (sector->xmax > a) {
                if (sector->xmin < a) // first
                    y = sector->propagate(E, y, sector->xmin - a, theta);

                if (sector->xmax >= b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, theta);
                    break;
                }

                y = sector->propagate(E, y, sector->xmax - sector->xmin, theta);
            }
        }
    } else {
        if (theta == 0)
            theta += M_PI;
        for (int i = sectorCount - 1; i >= 0; --i) {
            Sector *sector = sectors[i];
            if (sector->xmin < a) {
                if (sector->xmax > a) // first
                    y = sector->propagate(E, y, sector->xmin - a, theta);
                else
                    y = sector->propagate(E, y, sector->xmin - sector->xmax, theta);

                if (sector->xmin <= b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, theta);
                    break;
                }

            }
        }
    }
    return {y, theta};
}

tuple<double, double, double>
Matslise::calculateError(double E, const Y<> &left, const Y<> &right) const {
    Y<> l, r;
    double thetaL, thetaR;
    tie(l, thetaL) = propagate(E, left, xmin, match);
    tie(r, thetaR) = propagate(E, right, xmax, match);
    return make_tuple(l.y[1] * r.y[0] - r.y[1] * l.y[0],
                      l.dy[1] * r.y[0] + l.y[1] * r.dy[0] - (r.dy[1] * l.y[0] + r.y[1] * l.dy[0]),
                      thetaL - thetaR);
}

pair<int, double>
newtonIteration(const Matslise *ms, double E, const Y<> &left, const Y<> &right, double tol) {
    double adjust, error, derror, theta;
    int i = 0;
    do {
        tie(error, derror, theta) = ms->calculateError(E, left, right);
        adjust = error / derror;
        E = E - adjust;
        if (++i > 50) {
            cerr << "Newton-iteration did not converge for E=" << E << endl;
            break;
        }
    } while (fabs(adjust) > tol);

    int index = (int) round(theta / M_PI);
    if (index < 0)
        index = 0;
    return make_pair(index, E);
}

vector<pair<int, double>> *
Matslise::computeEigenvaluesByIndex(int Imin, int Imax, const Y<> &left,
                                    const Y<> &right) const {
    double Emin = -1;
    double Emax = 1;
    while (true) {
        double t = get<2>(calculateError(Emax, left, right)) / M_PI;
        int i = (int) floor(t);
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
            double t = get<2>(calculateError(Emin, left, right)) / M_PI;
            int i = (int) ceil(t);
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

vector<pair<int, double>> *
Matslise::computeEigenvalues(double Emin, double Emax, const Y<> &left, const Y<> &right) const {
    return computeEigenvalues(Emin, Emax, 0, INT_MAX, left, right);
};

vector<pair<int, double>> *
Matslise::computeEigenvalues(double Emin, double Emax, int Imin, int Imax, const Y<> &left,
                             const Y<> &right) const {
    if (Imin < 0)
        throw runtime_error("Matslise::computeEigenvalues(): Imin has to be at least 0");
    if (Imin > Imax)
        throw runtime_error("Matslise::computeEigenvalues(): Imax can't be less then Imin");
    vector<pair<int, double>> *eigenvalues = new vector<pair<int, double>>();
    queue<tuple<double, double, double, double, int>> toCheck;

    toCheck.push(make_tuple(Emin, get<2>(calculateError(Emin, left, right)) / M_PI,
                            Emax, get<2>(calculateError(Emax, left, right)) / M_PI,
                            0));

    double a, ta, b, tb, c, tc;
    int ia, ib, depth;
    while (!toCheck.empty()) {
        tie(a, ta, b, tb, depth) = toCheck.front();
        toCheck.pop();
        ia = (int) ceil(ta);
        ib = (int) ceil(tb);
        if (ta >= tb || ia == ib || ib <= Imin || Imax <= ia)
            continue;
        if (ia + 1 == ib)
            ++depth;

        c = (a + b) / 2;
        if (tb - ta < 0.05 || depth > 20)
            eigenvalues->push_back(newtonIteration(this, c, left, right, 1e-9));
        else {
            tc = get<2>(calculateError(c, left, right)) / M_PI;
            if (isnan(tc)) {
                cerr << "Matslise::computeEigenvalues(): some interval converted to NaN" << endl;
            } else {
                toCheck.push(make_tuple(a, ta, c, tc, depth));
                toCheck.push(make_tuple(c, tc, b, tb, depth));
            }
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

vector<Y<>> propagationSteps(const Matslise &ms, double E, const matslise::Y<> &left, const matslise::Y<> &right) {
    int n = ms.sectorCount;
    int m = n / 2 + 1;
    vector<Y<>> ys(n + 1);
    ys[0] = left;
    for (int i = 1; i <= m; ++i)
        ys[i] = ms.sectors[i - 1]->propagate(E, ys[i - 1], true);
    ys[n] = right;
    for (int i = n - 1; i > m; --i)
        ys[i] = ms.sectors[i]->propagate(E, ys[i + 1], false);
    Y<> yr = ms.sectors[m]->propagate(E, ys[m + 1], false);
    double s = ys[m].y[0] / yr.y[0];
    double norm = ys[m].dy[0] * ys[m].y[1] - ys[m].dy[1] * ys[m].y[0];
    norm -= s * s * (yr.dy[0] * yr.y[1] - yr.dy[1] * yr.y[0]);
    if (norm > 0) {
        norm = sqrt(norm);
    } else {
        cerr << "There are problems with the normalization." << endl;
        norm = 1;
    }
    int i = 0;
    for (; i <= m; ++i)
        ys[i] *= 1. / norm;
    for (; i <= n; ++i)
        ys[i] *= s / norm;
    return ys;
}

Array<Y<>, Dynamic, 1>
Matslise::computeEigenfunction(double E, const matslise::Y<> &left, const matslise::Y<> &right,
                               const ArrayXd &x) const {
    long n = x.size();
    for (int i = 1; i < n; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");
    if (x[0] < xmin || x[n - 1] > xmax)
        throw runtime_error("Matslise::computeEigenfunction(): x is out of range");

    vector<Y<>> steps = propagationSteps(*this, E, left, right);
    Array<Y<>, Dynamic, 1> ys(n);

    int sector = 0;
    for (int i = 0; i < n; ++i) {
        while (x[i] > sectors[sector]->xmax)
            ++sector;
        ys[i] = sectors[sector]->propagate(E, steps[sector], x[i] - sectors[sector]->xmin);
    }

    return ys;
}

std::function<Y<>(double)> Matslise::eigenfunctionCalculator(double E, const Y<> &left, const Y<> &right) const {
    vector<Y<>> ys = propagationSteps(*this, E, left, right);
    return [this, E, ys](double x) -> Y<> {
        int a = 0;
        int b = this->sectorCount;
        while (a + 1 < b) {
            int c = (a + b) / 2;
            if (x < this->sectors[c]->xmin)
                b = c;
            else
                a = c;
        }
        return this->sectors[a]->propagate(E, ys[a], x - this->sectors[a]->xmin);
    };
}