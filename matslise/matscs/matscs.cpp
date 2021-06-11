#include <functional>
#include "../matscs.h"
#include "../util/find_sector.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;
using namespace Eigen;

template<typename Scalar>
template<int r>
Y<Scalar, Dynamic, r>
Matscs<Scalar>::propagateColumn(const Scalar &E, const Y<Scalar, Dynamic, r> &_y,
                                const Scalar &a, const Scalar &b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Scalar, Dynamic, r> y = _y;
    int sectorIndex = find_sector<Matscs<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;

    do {
        const value_ptr<Sector> &sector = sectors[sectorIndex];
        y = sector->propagateColumn(E, y, a, b, use_h);
        if (sector->contains(b))
            break;
        sectorIndex += direction;
    } while (sectorIndex >= 0 && sectorIndex < (long) sectors.size());

    return y;
}

template<typename Scalar>
pair<Y<Scalar, Dynamic>, Scalar>
Matscs<Scalar>::propagate(const Scalar &E, const Y<Scalar, Dynamic> &_y,
                          const Scalar &a, const Scalar &b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Scalar, Dynamic> y = _y;
    int sectorIndex = find_sector<Matscs<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Scalar dTheta;
    Scalar argdet = 0;
    do {
        const value_ptr<Sector> &sector = sectors[sectorIndex];
        tie(y, dTheta) = sector->propagate(E, y, a, b, use_h);
        argdet += dTheta;
        if (sector->contains(b))
            break;
        sectorIndex += direction;
    } while (sectorIndex >= 0 && sectorIndex < (long) sectors.size());

    return {y, argdet};
}

template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic> Matscs<Scalar>::propagatePsi(
        const Scalar &E, const Matrix<Scalar, Dynamic, Dynamic> &_psi, const Scalar &a, const Scalar &b) const {
    Matrix<Scalar, Dynamic, Dynamic> psi = _psi;

    int sectorIndex = find_sector<Matscs<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;

    do {
        const value_ptr<Sector> &sector = sectors[sectorIndex];
        psi = sector->propagatePsi(E, psi, a, b);
        if (sector->contains(b))
            break;
        sectorIndex += direction;
    } while (sectorIndex >= 0 && sectorIndex < (long) sectors.size());

    return psi;
}

template<typename Scalar>
vector<Y<Scalar, Dynamic>> *Matscs<Scalar>::eigenfunction(const Scalar &E, vector<Scalar> &x) const {
    sort(x.begin(), x.end());
    auto *ys = new vector<Y<Scalar, Dynamic>>();

    auto iterator = x.begin();

    while (iterator != x.end() && *iterator < xmin - EPS)
        iterator = x.erase(iterator);

    Y<Scalar, Dynamic> y(n);
    for (int i = 0; iterator != x.end(); ++iterator) {
        const value_ptr<Sector> &sector = sectors[i];
        while (*iterator > sector->max) {
            y = sector->calculateT(E) * y;
            ++i;
            if (i >= sectorCount)
                goto allSectorsDone;
        }

        ys->push_back(sector->calculateT(E, *iterator - sector->min) * y);
    }

    allSectorsDone:
    while (iterator != x.end() && *iterator > xmax + EPS)
        iterator = x.erase(iterator);

    return ys;
}


#include "../util/sectorbuilder.impl.h"

#define INSTANTIATE_PROPAGATE(Scalar, r)\
template Y<Scalar, Dynamic, r>\
Matscs<Scalar>::propagateColumn<r>(\
    const Scalar &E, const Y<Scalar, Dynamic, r> &y, const Scalar &a, const Scalar &b, bool use_h) const;

#define INSTANTIATE_MORE(Scalar)\
INSTANTIATE_PROPAGATE(Scalar, 1)\
INSTANTIATE_PROPAGATE(Scalar, -1)\
INSTANTIATE_SECTOR_BUILDER(Matscs<Scalar>)

#include "instantiate.h"
