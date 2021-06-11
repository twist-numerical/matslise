#include <array>
#include <vector>
#include <queue>
#include "../matslise.h"
#include "matslise_formulas.h"
#include "../util/theta.h"
#include "../util/constants.h"
#include "../util/find_sector.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
pair<Y<Scalar>, Scalar>
Matslise<Scalar>::propagate(const Scalar &E, const Y<Scalar> &_y, const Scalar &a, const Scalar &b, bool use_h) const {
    if (!domain.contains(a) || !domain.contains(b))
        throw runtime_error("Matslise::propagate(): a and b should be in the interval");
    Y<Scalar> y = _y;
    Scalar dTheta;
    long sectorIndex = find_sector<Matslise<Scalar>>(this, a);
    Scalar theta = sectors[sectorIndex]->theta0(E, y);
    if (theta < 0 || (a > b && theta == 0))
        theta += constants<Scalar>::PI;
    do {
        const value_ptr<Sector> &sector = sectors[sectorIndex];
        tie(y, dTheta) = sector->template propagate<true>(E, y, a, b, use_h);
        theta += dTheta;
        if (sector->contains(b))
            break;
        sectorIndex += a < b ? 1 : -1;
    } while (sectorIndex >= 0 && sectorIndex < (long) sectors.size());
    return {y, theta};
}

template<typename Scalar>
tuple<Scalar, Scalar, Scalar>
Matslise<Scalar>::matchingError(const Scalar &E, const Y<Scalar> &left, const Y<Scalar> &right, bool use_h) const {
    Y<Scalar> l, r;
    Scalar thetaL, thetaR;
    tie(l, thetaL) = propagate(E, left, domain.min(), sectors[matchIndex]->max, use_h);
    tie(r, thetaR) = propagate(E, right, domain.max(), sectors[matchIndex]->max, use_h);
    return make_tuple(l.data[1] * r.data[0] - r.data[1] * l.data[0],
                      l.data[3] * r.data[0] + l.data[1] * r.data[2] - (r.data[3] * l.data[0] + r.data[1] * l.data[2]),
                      thetaL - thetaR);
}

template<typename Scalar>
pair<int, Scalar>
newtonIteration(const Matslise<Scalar> *ms, Scalar E, const Y<Scalar> &left, const Y<Scalar> &right, bool use_h) {
    Scalar adjust, error, derror, theta;
    int i = 0;
    do {
        tie(error, derror, theta) = ms->matchingError(E, left, right, use_h);
        adjust = error / derror;
        E -= adjust;
    } while (++i < 20 && abs(adjust) > ms->tolerance);

    if (i >= 20 && abs(adjust) > 100 * ms->tolerance) {
        // cerr << "Newton-iteration did not converge for E=" << (double) E << endl;
    }

    int index = (int) round(theta / constants<Scalar>::PI);
    if (index < 0)
        index = 0;
    return make_pair(index, E);
}

template<typename Scalar>
Scalar Matslise<Scalar>::estimatePotentialMinimum() const {
    auto iterator = this->sectors.begin();
    Scalar minimal = (*iterator++)->vs[0];
    for (; iterator != this->sectors.end(); ++iterator)
        minimal = min(minimal, (*iterator)->vs[0]);
    return minimal;
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

    toCheck.push(make_tuple(Emin, get<2>(ms->matchingError(Emin, left, right)) / constants<Scalar>::PI,
                            Emax, get<2>(ms->matchingError(Emax, left, right)) / constants<Scalar>::PI,
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

        c = ia + 1 < ib || tb - ta < 1e-5 || depth % 2 == 0
            ? .5 * (a + b)
            : ((tb - ia) * a - (ta - ia) * b) / (tb - ta);
        if ((tb - ta < 0.01 && depth > 3) || depth > 30) {
            eigenvalues.push_back(newtonIteration<Scalar>(ms, c, left, right, true));
        } else {
            tc = get<2>(ms->matchingError(c, left, right)) / constants<Scalar>::PI;
            if (isnan(tc)) {
                // cerr << "Matslise::computeEigenvalues(): some interval converted to NaN" << endl;
            } else if (ia + 1 < ib) {
                toCheck.push(make_tuple(a, ta, c, tc, 1 + depth));
                toCheck.push(make_tuple(c, tc, b, tb, 1 + depth));
            } else {
                if (abs(tc - ia) < 1e-8)
                    eigenvalues.push_back(newtonIteration<Scalar>(ms, c, left, right, true));
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
Matslise<Scalar>::eigenvaluesByIndex(int Imin, int Imax, const Y<Scalar> &left, const Y<Scalar> &right) const {
    Scalar Emin = -1;
    Scalar Emax = 1;
    while (true) {
        Scalar t = get<2>(matchingError(Emax, left, right)) / constants<Scalar>::PI;
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
            Scalar t = get<2>(matchingError(Emin, left, right)) / constants<Scalar>::PI;
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
Matslise<Scalar>::eigenvalues(const Scalar &Emin, const Scalar &Emax, const Y<Scalar> &left,
                              const Y<Scalar> &right) const {
    return computeEigenvaluesHelper(this, Emin, Emax, 0, INT_MAX, left, right);
}

template<typename Scalar>
Scalar
Matslise<Scalar>::eigenvalueError(const Scalar &E, const Y<Scalar> &left, const Y<Scalar> &right, int) const {
    return abs(E - newtonIteration<Scalar>(this, E, left, right, false).second);
}

template<typename Scalar>
vector<Y<Scalar>> propagationSteps(const Matslise<Scalar> &ms, Scalar E,
                                   const Y<Scalar> &left, const Y<Scalar> &right) {
    vector<Y<Scalar>> ys(ms.sectorCount + 1);
    ys[0] = left;
    for (int i = 0; i <= ms.matchIndex; ++i) {
        const value_ptr<typename Matslise<Scalar>::Sector> &sector = ms.sectors[i];
        ys[i + 1] = sector->template propagate<false>(E, ys[i], sector->min, sector->max);
    }
    Y<Scalar> yl = ys[ms.matchIndex + 1];

    ys[ms.sectorCount] = right;
    for (int i = ms.sectorCount - 1; i > ms.matchIndex; --i) {
        const value_ptr<typename Matslise<Scalar>::Sector> &sector = ms.sectors[i];
        ys[i] = sector->template propagate<false>(E, ys[i + 1], sector->max, sector->min);
    }
    Y<Scalar> &yr = ys[ms.matchIndex + 1];

    Scalar s = (abs(yr.data[0]) + abs(yl.data[0]) > abs(yr.data[1]) + abs(yl.data[1])) ? yl.data[0] / yr.data[0] :
               yl.data[1] / yr.data[1];
    Scalar norm = yl.data[2] * yl.data[1] - yl.data[3] * yl.data[0];
    norm -= s * s * (yr.data[2] * yr.data[1] - yr.data[3] * yr.data[0]);
    if (norm > 0) {
        norm = sqrt(norm);
    } else {
        cerr << "There are problems with the normalization." << endl;
        norm = 1;
    }
    int i = 0;
    for (; i <= ms.matchIndex; ++i)
        ys[i] *= ((Scalar) 1.) / norm;
    for (; i <= ms.sectorCount; ++i)
        ys[i] *= s / norm;
    return ys;
}

template<typename Scalar, typename Iterator>
Iterator findSector(const Scalar &x, Iterator a, Iterator b) {
    while (b - a > 1) {
        Iterator c = a + (b - a) / 2;
        if (x < (**c).min)
            b = c;
        else
            a = c;
    }
    return a;
}

template<typename Scalar>
typename Matslise<Scalar>::Eigenfunction Matslise<Scalar>::eigenfunction(
        const Scalar &E, const Y<Scalar> &left, const Y<Scalar> &right, int) const {
    shared_ptr<vector<Y<Scalar >>> steps = make_shared<vector<Y<Scalar >>>(
            std::move(propagationSteps(*this, E, left, right)));

    return {[this, E, steps](const Scalar &x) {
        if (x < domain.template min<0>() || x > domain.template max<0>())
            return Y<Scalar>();
        auto sector = findSector(x, sectors.begin(), sectors.end());
        int sectorIndex = sector - sectors.begin();

        if ((**sector).direction == forward) {
            return (**sector).template propagate<false>(E, (*steps)[sectorIndex], (**sector).min, x);
        } else {
            return (**sector).template propagate<false>(E, (*steps)[sectorIndex + 1], (**sector).max, x);
        }
    }, [this, E, steps](const Array<Scalar, Dynamic, 1> &x) {
        Eigen::Index n = x.size();

        Array<Y<Scalar>, Dynamic, 1> ys(n);

        int sectorIndex = 0;
        auto sector = sectors.begin();
        for (Eigen::Index i = 0; i < n; ++i) {
            if (x[i] < domain.template min<0>() || x[i] > domain.template max<0>()) {
                ys[i] = Y<Scalar>();
                continue;
            }
            if (x[i] < (**sector).min) {
                sector = findSector(x[i], sectors.begin(), sector);
                sectorIndex = sector - sectors.begin();
            } else if (x[i] > (**sector).max) {
                sector = findSector(x[i], sector + 1, sectors.end());
                sectorIndex = sector - sectors.begin();
            }
            if ((**sector).direction == forward) {
                ys[i] = (**sector).template propagate<false>(E, (*steps)[sectorIndex], (**sector).min, x[i]);
            } else {
                ys[i] = (**sector).template propagate<false>(E, (*steps)[sectorIndex + 1], (**sector).max, x[i]);
            }

        }
        return ys;
    }};
}

#include "../util/instantiate.h"