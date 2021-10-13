#include "../matslise.h"
#include "../util/find_sector.h"

using namespace matslise;
using namespace Eigen;
using namespace std;

template<typename Scalar>
Y<Scalar, 1, 2>
PeriodicMatslise<Scalar>::propagate(
        const Scalar &E, const Y<Scalar, 1, 2> &_y,
        const Scalar &a, const Scalar &b, bool use_h) const {
    Y<Scalar, 1, 2> y = _y;

    long sectorIndex = find_sector<Matslise<Scalar>>(&matslise, a);
    do {
        const value_ptr<typename Matslise<Scalar>::Sector> &sector = matslise.sectors[sectorIndex];
        y = sector->template propagate<false, 2>(E, y, a, b, use_h);
        if (sector->contains(b))
            break;
        sectorIndex += a < b ? 1 : -1;
    } while (sectorIndex >= 0 && sectorIndex < (long) matslise.sectors.size());
    return y;
}

template<typename Scalar>
pair<Scalar, Scalar> PeriodicMatslise<Scalar>::matchingError(const Scalar &E, bool use_h) const {
    Y<Scalar, 1, 2> l = Y<Scalar, 1, 2>::Periodic();
    Y<Scalar, 1, 2> r = Y<Scalar, 1, 2>::Periodic();
    l = propagate(E, l, matslise.domain.min(), matslise.sectors[matslise.matchIndex]->max, use_h);
    r = propagate(E, r, matslise.domain.max(), matslise.sectors[matslise.matchIndex]->max, use_h);
    Matrix<Scalar, 2, 2> error = l.y() - r.y();
    Matrix<Scalar, 2, 2> dError = l.ydE() - r.ydE();
    return make_pair(error.determinant(), 0);
}

#include "instantiate.h"
