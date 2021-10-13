#include "../matslise.h"
#include "../util/find_sector.h"

using namespace matslise;
using namespace Eigen;
using namespace std;

template<typename Scalar>
pair<Y<Scalar, 1, 2>, Array<Scalar, 2, 1>>
PeriodicMatslise<Scalar>::propagate(
        const Scalar &E, const Y<Scalar, 1, 2> &y,
        const Scalar &a, const Scalar &b, bool use_h) const {
    return matslise.template propagate<2>(E, y, a, b, use_h);
}

template<typename Scalar>
tuple<Matrix<Scalar, 2, 2>, Matrix<Scalar, 2, 2>, Array<Scalar, 2, 1>>
PeriodicMatslise<Scalar>::matchingError(const Scalar &E, bool use_h) const {
    Y<Scalar, 1, 2> l = Y<Scalar, 1, 2>::Periodic();
    Y<Scalar, 1, 2> r = Y<Scalar, 1, 2>::Periodic();
    Array<Scalar, 2, 1> thetaL, thetaR;
    tie(l, thetaL) = propagate(E, l, matslise.domain.min(), matslise.sectors[matslise.matchIndex]->max, use_h);
    tie(r, thetaR) = propagate(E, r, matslise.domain.max(), matslise.sectors[matslise.matchIndex]->max, use_h);
    Matrix<Scalar, 2, 2> error = l.y() - r.y();
    Matrix<Scalar, 2, 2> dError = l.ydE() - r.ydE();
    return make_tuple(error, dError, thetaL - thetaR);
}

#include "instantiate.h"
