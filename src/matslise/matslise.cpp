//
// Created by toon on 5/16/18.
//

#include <array>
#include <vector>
#include <queue>
#include "../matslise.h"
#include "matslise_formulas.h"
#include "../util/theta.h"
#include "../util/fmath.h"
#include "../util/find_sector.h"

#define EPS (1.e-12)

using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
pair<Y<Scalar>, Scalar>
Matslise<Scalar>::propagate(Scalar E, const Y<Scalar> &_y, Scalar a, Scalar b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matslise::propagate(): a and b should be in the interval");
    Y<Scalar> y = _y;
    Scalar theta = matslise::theta(y);
    if (a == xmax && theta == 0)
        theta += M_PI;
    int sectorIndex = find_sector<Matslise<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        y = sector->propagate(E, y, a, b, theta, use_h);
        sectorIndex += direction;
    } while (!sector->contains(b));
    return {y, theta};
}

template<typename Scalar>
tuple<Scalar, Scalar, Scalar>
Matslise<Scalar>::calculateError(Scalar E, const Y<Scalar> &left, const Y<Scalar> &right, bool use_h) const {
    Y<Scalar> l, r;
    Scalar thetaL, thetaR;
    tie(l, thetaL) = propagate(E, left, xmin, match, use_h);
    tie(r, thetaR) = propagate(E, right, xmax, match, use_h);
    return make_tuple(l.y[1] * r.y[0] - r.y[1] * l.y[0],
                      l.dy[1] * r.y[0] + l.y[1] * r.dy[0] - (r.dy[1] * l.y[0] + r.y[1] * l.dy[0]),
                      thetaL - thetaR);
}

template<typename Scalar>
pair<int, Scalar>
newtonIteration(const Matslise<Scalar> *ms, Scalar E, const Y<Scalar> &left, const Y<Scalar> &right, Scalar tol,
                bool use_h) {
    Scalar adjust, error, derror, theta;
    int i = 0;
    do {
        tie(error, derror, theta) = ms->calculateError(E, left, right, use_h);
        adjust = error / derror;
        E = E - adjust;
        if (++i > 50) {
            cerr << "Newton-iteration did not converge for E=" << (double) E << endl;
            break;
        }
    } while (fmath<Scalar>::abs(adjust) > tol);

    int index = fmath<Scalar>::round(theta / M_PI);
    if (index < 0)
        index = 0;
    return make_pair(index, E);
}

template<typename Scalar>
vector<pair<int, Scalar>> *
computeEigenvaluesHelper(const Matslise<Scalar> *ms, Scalar Emin, Scalar Emax, int Imin, int Imax,
                         const Y<Scalar> &left, const Y<Scalar> &right) {
    if (Imin < 0)
        throw runtime_error("Matslise::computeEigenvalues(): Imin has to be at least 0");
    if (Imin > Imax)
        throw runtime_error("Matslise::computeEigenvalues(): Imax can't be less then Imin");
    vector<pair<int, Scalar>> *eigenvalues = new vector<pair<int, Scalar>>();
    queue<tuple<Scalar, Scalar, Scalar, Scalar, int>> toCheck;

    toCheck.push(make_tuple(Emin, get<2>(ms->calculateError(Emin, left, right)) / M_PI,
                            Emax, get<2>(ms->calculateError(Emax, left, right)) / M_PI,
                            0));

    Scalar a, ta, b, tb, c, tc;
    int ia, ib, depth;
    while (!toCheck.empty()) {
        tie(a, ta, b, tb, depth) = toCheck.front();
        toCheck.pop();
        ia = fmath<Scalar>::ceil(ta);
        ib = fmath<Scalar>::ceil(tb);
        if (ta >= tb || ia == ib || ib <= Imin || Imax <= ia)
            continue;
        if (ia + 1 == ib)
            ++depth;

        c = (a + b) / 2;
        if (tb - ta < 0.05 || depth > 20)
            eigenvalues->push_back(newtonIteration<Scalar>(ms, c, left, right, 1e-9, true));
        else {
            tc = get<2>(ms->calculateError(c, left, right)) / M_PI;
            if (fmath<Scalar>::isnan(tc)) {
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


template<typename Scalar>
vector<pair<int, Scalar>> *
Matslise<Scalar>::computeEigenvaluesByIndex(int Imin, int Imax, const Y<Scalar> &left, const Y<Scalar> &right) const {
    Scalar Emin = -1;
    Scalar Emax = 1;
    while (true) {
        Scalar t = get<2>(calculateError(Emax, left, right)) / M_PI;
        int i = fmath<Scalar>::floor(t);
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
            Scalar t = get<2>(calculateError(Emin, left, right)) / M_PI;
            int i = fmath<Scalar>::ceil(t);
            if (i <= Imin)
                break;
            else {
                if (i > Imax)
                    Emax = Emin;
                Emin *= 2;
            }
        }
    }
    return computeEigenvaluesHelper(this, Emin, Emax, Imin, Imax, left, right);
}

template<typename Scalar>
vector<pair<int, Scalar>> *
Matslise<Scalar>::computeEigenvalues(Scalar Emin, Scalar Emax, const Y<Scalar> &left, const Y<Scalar> &right) const {
    return computeEigenvaluesHelper(this, Emin, Emax, 0, INT_MAX, left, right);
}

template<typename Scalar>
Scalar Matslise<Scalar>::computeEigenvalueError(Scalar E, const Y<Scalar> &left, const Y<Scalar> &right) const {
    return fmath<Scalar>::abs(E - newtonIteration<Scalar>(this, E, left, right, 1e-9, false).second);
}

template<typename Scalar>
Matslise<Scalar>::~Matslise() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
}

template<typename Scalar>
vector<Y<Scalar>> propagationSteps(const Matslise<Scalar> &ms, Scalar E,
                                   const Y<Scalar> &left, const Y<Scalar> &right) {
    unsigned int n = static_cast<unsigned int>(ms.sectorCount);
    vector<Y<Scalar>> ys(n + 1);
    ys[0] = left;
    unsigned int m = 0;
    for (unsigned int i = 1; ms.sectors[i - 1]->min < ms.match; ++i) {
        ys[i] = ms.sectors[i - 1]->propagate(E, ys[i - 1], true);
        m = i;
    }
    ys[n] = right;
    for (unsigned int i = n - 1; i > m; --i)
        ys[i] = ms.sectors[i]->propagate(E, ys[i + 1], false);
    Y<Scalar> yr = ms.sectors[m]->propagate(E, ys[m + 1], false);
    Y<Scalar> &yl = ys[m];
    Scalar s = (fmath<Scalar>::abs(yr.y[0]) + fmath<Scalar>::abs(yl.y[0]) >
                fmath<Scalar>::abs(yr.y[1]) + fmath<Scalar>::abs(yl.y[1])) ? yl.y[0] / yr.y[0] : yl.y[1] / yr.y[1];
    Scalar norm = yl.dy[0] * yl.y[1] - yl.dy[1] * yl.y[0];
    norm -= s * s * (yr.dy[0] * yr.y[1] - yr.dy[1] * yr.y[0]);
    if (norm > 0) {
        norm = fmath<Scalar>::sqrt(norm);
    } else {
        cerr << "There are problems with the normalization." << endl;
        norm = 1;
    }
    unsigned int i = 0;
    for (; i <= m; ++i)
        ys[i] *= ((Scalar) 1.) / norm;
    for (; i <= n; ++i)
        ys[i] *= s / norm;
    return ys;
}

template<typename Scalar>
Array<Y<Scalar>, Dynamic, 1>
Matslise<Scalar>::computeEigenfunction(Scalar E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                                       const Array<Scalar, Dynamic, 1> &x) const {
    long n = x.size();
    for (int i = 1; i < n; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");
    if (x[0] < xmin || x[n - 1] > xmax)
        throw runtime_error("Matslise::computeEigenfunction(): x is out of range");

    vector<Y<Scalar>> steps = propagationSteps(*this, E, left, right);
    Array<Y<Scalar>, Dynamic, 1> ys(n);

    unsigned int sector = 0;
    Scalar theta = 0;
    for (int i = 0; i < n; ++i) {
        while (x[i] > sectors[sector]->max)
            ++sector;
        if (sectors[sector]->backward) {
            ys[i] = sectors[sector]->propagate(E, steps[sector + 1], sectors[sector]->max, x[i], theta);
        } else {
            ys[i] = sectors[sector]->propagate(E, steps[sector], sectors[sector]->min, x[i], theta);
        }

    }

    return ys;
}

template<typename Scalar>
std::function<Y<Scalar>(Scalar)> Matslise<Scalar>::eigenfunctionCalculator(
        Scalar E, const Y<Scalar> &left, const Y<Scalar> &right) const {
    vector<Y<Scalar>> ys = propagationSteps(*this, E, left, right);
    return [this, E, ys](Scalar x) -> Y<Scalar> {
        int a = 0;
        int b = this->sectorCount;
        while (a + 1 < b) {
            int c = (a + b) / 2;
            if (x < this->sectors[c]->min)
                b = c;
            else
                a = c;
        }
        const Matslise<Scalar>::Sector *sector = this->sectors[a];
        Scalar theta = 0;
        if (sector->backward) {
            return sector->propagate(E, ys[static_cast<unsigned long>(a + 1)], sector->max, x, theta);
        } else {
            return sector->propagate(E, ys[static_cast<unsigned long>(a)], sector->min, x, theta);
        }
    };
}

#include "../util/instantiate.h"