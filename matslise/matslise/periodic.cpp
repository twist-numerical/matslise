#include "../matslise.h"
#include "../util/find_sector.h"
#include "../util/constants.h"

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
tuple<Scalar, Scalar, Array<Scalar, 2, 1>>
PeriodicMatslise<Scalar>::matchingError(const Scalar &E, bool use_h) const {
    Y<Scalar, 1, 2> l = Y<Scalar, 1, 2>::Periodic();
    Y<Scalar, 1, 2> r = Y<Scalar, 1, 2>::Periodic();
    Array<Scalar, 2, 1> thetaL, thetaR;
    tie(l, thetaL) = propagate(E, l, matslise.domain.min(), matslise.sectors[matslise.matchIndex]->max, use_h);
    tie(r, thetaR) = propagate(E, r, matslise.domain.max(), matslise.sectors[matslise.matchIndex]->max, use_h);
    Matrix<Scalar, 2, 2> err = l.y() - r.y();
    Matrix<Scalar, 2, 2> dErr = l.ydE() - r.ydE();
    Scalar error = err(0, 0) * err(1, 1) - err(0, 1) * err(1, 0);
    Scalar dError = dErr(0, 0) * err(1, 1) + err(0, 0) * dErr(1, 1) - dErr(0, 1) * err(1, 0) - err(0, 1) * dErr(1, 0);

    Array<Scalar, 2, 1> theta = thetaL - thetaR;
    theta /= constants<Scalar>::PI;
    theta[1] += 1; // adjust for initial conditions
    return make_tuple(error, dError, theta);
}

template<typename Scalar>
vector<unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>
PeriodicMatslise<Scalar>::eigenfunction(const Scalar &E) const {
    Array<Scalar, 2, 1> thetaL, thetaR;
    Scalar &match = matslise.sectors[matslise.matchIndex]->max;
    Y<Scalar, 1, 2> l = propagate(E, Y<Scalar, 1, 2>::Periodic(), matslise.domain.min(), match).first;
    Y<Scalar, 1, 2> r = propagate(E, Y<Scalar, 1, 2>::Periodic(), matslise.domain.max(), match).first;
    Matrix<Scalar, 2, 2> err = (l.y() - r.y());
    cout << err << endl;

    EigenSolver<Matrix<Scalar, 2, 2>> solver(err, true);
    Matrix<complex<Scalar>, 2, 1> eigenvalues = solver.eigenvalues();
    Matrix<complex<Scalar>, 2, 2> eigenvectors = solver.eigenvectors();

    vector<unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>> result;
    cout << "eigenvalues: "<< eigenvalues.transpose() << endl;
    cout << eigenvectors << endl;
    cout << endl;
    result.reserve(2);
    if (norm(eigenvalues[0]) < 1e-4) {
        Y<> boundary;
        boundary.y() = eigenvectors.col(0).real();
        result.emplace_back(matslise.eigenfunction(E, boundary, boundary));
    }
    if (norm(eigenvalues[1]) < 1e-4) {
        Y<> boundary;
        boundary.y() = eigenvectors.col(1).real();
        result.emplace_back(matslise.eigenfunction(E, boundary, boundary));
    }

    return result;
}

#include "instantiate.h"
