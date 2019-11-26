#include <array>
#include <vector>
#include <queue>
#include "../matslise.h"
#include "matslise_formulas.h"
#include "../util/theta.h"
#include "../util/constants.h"
#include "../util/find_sector.h"
#include "../util/addTheta.h"

#define EPS (1.e-12)

using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
pair<Y<Scalar>, Scalar>
Matslise<Scalar>::propagate(const Scalar &E, const Y<Scalar> &_y, const Scalar &a, const Scalar &b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matslise::propagate(): a and b should be in the interval");
    Y<Scalar> y = _y;
    Scalar theta = matslise::theta(y);
    if (a == xmax && theta == 0)
        theta += constants<Scalar>::PI;
    int sectorIndex = find_sector<Matslise<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        y = addTheta(sector->propagate(E, y, a, b, use_h), theta);
        sectorIndex += direction;
    } while (!sector->contains(b));
    return {y, theta};
}

template<typename Scalar>
tuple<Scalar, Scalar, Scalar>
Matslise<Scalar>::calculateError(const Scalar &E, const Y<Scalar> &left, const Y<Scalar> &right, bool use_h) const {
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
        E -= adjust;
    } while (++i < 50 && abs(adjust) > tol);

    if (i >= 50) {
        cerr << "Newton-iteration did not converge for E=" << (double) E << endl;
    }

    int index = (int) round(theta / constants<Scalar>::PI);
    if (index < 0)
        index = 0;
    return make_pair(index, E);
}

template<typename Scalar>
vector<pair<int, Scalar>>
computeEigenvaluesHelper(const Matslise<Scalar> *ms, Scalar Emin, Scalar Emax, int Imin, int Imax,
                         const Y<Scalar> &left, const Y<Scalar> &right) {
    if (Imin < 0)
        throw runtime_error("Matslise::computeEigenvalues(): Imin has to be at least 0");
    if (Imin > Imax)
        throw runtime_error("Matslise::computeEigenvalues(): Imax can't be less then Imin");
    vector<pair<int, Scalar>> eigenvalues;
    queue<tuple<Scalar, Scalar, Scalar, Scalar, int>> toCheck;

    toCheck.push(make_tuple(Emin, get<2>(ms->calculateError(Emin, left, right)) / constants<Scalar>::PI,
                            Emax, get<2>(ms->calculateError(Emax, left, right)) / constants<Scalar>::PI,
                            0));

    Scalar a, ta, b, tb, c, tc;
    int ia, ib, depth;
    while (!toCheck.empty()) {
        tie(a, ta, b, tb, depth) = toCheck.front();
        toCheck.pop();
        ia = (int) ceil(ta);
        ib = (int) ceil(tb);
        if (ta >= tb || ia == ib || ib <= Imin || Imax <= ia)
            continue;

        c = ia + 1 < ib || depth % 2 == 0 ? .5 * (a + b) : ((tb - ia) * a - (ta - ia) * b) / (tb - ta);
        if ((tb - ta < 0.05 && depth > 3) || depth > 20) {
            eigenvalues.push_back(newtonIteration<Scalar>(ms, c, left, right, 1e-9, true));
        } else {
            tc = get<2>(ms->calculateError(c, left, right)) / constants<Scalar>::PI;
            if (isnan(tc)) {
                cerr << "Matslise::computeEigenvalues(): some interval converted to NaN" << endl;
            } else if (ia + 1 < ib) {
                toCheck.push(make_tuple(a, ta, c, tc, depth));
                toCheck.push(make_tuple(c, tc, b, tb, depth));
            } else {
                if (abs(tc - ia) < 1e-8)
                    eigenvalues.push_back(newtonIteration<Scalar>(ms, c, left, right, 1e-9, true));
                else if ((ta - ia) * (tc - ia) < 0)
                    toCheck.push(make_tuple(a, ta, c, tc, depth + 1));
                else
                    toCheck.push(make_tuple(c, tc, b, tb, depth + 1));
            }
        }
    }

    sort(eigenvalues.begin(), eigenvalues.end());

    return eigenvalues;
}


template<typename Scalar>
vector<pair<int, Scalar>>
Matslise<Scalar>::computeEigenvaluesByIndex(int Imin, int Imax, const Y<Scalar> &left, const Y<Scalar> &right) const {
    Scalar Emin = -1;
    Scalar Emax = 1;
    while (true) {
        Scalar t = get<2>(calculateError(Emax, left, right)) / constants<Scalar>::PI;
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
            Scalar t = get<2>(calculateError(Emin, left, right)) / constants<Scalar>::PI;
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
    return computeEigenvaluesHelper(this, Emin, Emax, Imin, Imax, left, right);
}

template<typename Scalar>
vector<pair<int, Scalar>>
Matslise<Scalar>::computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const Y<Scalar> &left,
                                     const Y<Scalar> &right) const {
    return computeEigenvaluesHelper(this, Emin, Emax, 0, INT_MAX, left, right);
}

template<typename Scalar>
Scalar
Matslise<Scalar>::computeEigenvalueError(const Scalar &E, const Y<Scalar> &left, const Y<Scalar> &right, int) const {
    return abs(E - newtonIteration<Scalar>(this, E, left, right, 1e-9, false).second);
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
    for (unsigned int i = 1; ms.match > ms.sectors[i - 1]->min; ++i) {
        ys[i] = ms.sectors[i - 1]->propagateForward(E, ys[i - 1]).first;
        m = i;
    }
    ys[n] = right;
    for (unsigned int i = n - 1; i > m; --i)
        ys[i] = ms.sectors[i]->propagateBackward(E, ys[i + 1]).first;
    Y<Scalar> yr = ms.sectors[m]->propagateBackward(E, ys[m + 1]).first;
    Y<Scalar> &yl = ys[m];
    Scalar s = (abs(yr.y[0]) + abs(yl.y[0]) > abs(yr.y[1]) + abs(yl.y[1])) ? yl.y[0] / yr.y[0] : yl.y[1] / yr.y[1];
    Scalar norm = yl.dy[0] * yl.y[1] - yl.dy[1] * yl.y[0];
    norm -= s * s * (yr.dy[0] * yr.y[1] - yr.dy[1] * yr.y[0]);
    if (norm > 0) {
        norm = sqrt(norm);
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
Matslise<Scalar>::computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &left,
                                       const matslise::Y<Scalar> &right,
                                       const Array<Scalar, Dynamic, 1> &x, int) const {
    Eigen::Index n = x.size();
    for (int i = 1; i < n; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");
    if (x[0] < xmin || x[n - 1] > xmax)
        throw runtime_error("Matslise::computeEigenfunction(): x is out of range");

    vector<Y<Scalar>> steps = propagationSteps(*this, E, left, right);
    Array<Y<Scalar>, Dynamic, 1> ys(n);

    unsigned int sector = 0;
    for (int i = 0; i < n; ++i) {
        while (x[i] > sectors[sector]->max)
            ++sector;
        if (sectors[sector]->backward) {
            ys[i] = sectors[sector]->propagate(E, steps[sector + 1], sectors[sector]->max, x[i]).first;
        } else {
            ys[i] = sectors[sector]->propagate(E, steps[sector], sectors[sector]->min, x[i]).first;
        }

    }

    return ys;
}

template<typename Scalar>
std::function<Y<Scalar>(Scalar)> Matslise<Scalar>::eigenfunctionCalculator(
        const Scalar &E, const Y<Scalar> &left, const Y<Scalar> &right, int) const {
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
        if (sector->backward) {
            return sector->propagate(E, ys[static_cast<unsigned long>(a + 1)], sector->max, x).first;
        } else {
            return sector->propagate(E, ys[static_cast<unsigned long>(a)], sector->min, x).first;
        }
    };
}

#include "../util/instantiate.h"