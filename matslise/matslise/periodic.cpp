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
pair<Y<Scalar, 1, 2>, Array<Scalar, 2, 1>>
PeriodicMatslise<Scalar>::matchingError(const Scalar &E, bool use_h) const {
    Y<Scalar, 1, 2> l = Y<Scalar, 1, 2>::Periodic();
    Y<Scalar, 1, 2> r = Y<Scalar, 1, 2>::Periodic();
    Array<Scalar, 2, 1> thetaL, thetaR;
    tie(l, thetaL) = propagate(E, l, matslise.domain.min(), matslise.sectors[matslise.matchIndex]->max, use_h);
    tie(r, thetaR) = propagate(E, r, matslise.domain.max(), matslise.sectors[matslise.matchIndex]->max, use_h);

    Y<Scalar, 1, 2> error = l;
    error.data -= r.data;
    Array<Scalar, 2, 1> theta = thetaL - thetaR;
    theta /= constants<Scalar>::PI;
    theta[1] += 1; // adjust for initial conditions
    return {error, theta};
}

template<typename Scalar>
pair<Scalar, Scalar> calculateMatchingError(const Y<Scalar, 1, 2> &yError) {
    auto err = yError.y();
    auto dErr = yError.ydE();
    Scalar error = err(0, 0) * err(1, 1) - err(0, 1) * err(1, 0);
    Scalar dError = dErr(0, 0) * err(1, 1) + err(0, 0) * dErr(1, 1) - dErr(0, 1) * err(1, 0) - err(0, 1) * dErr(1, 0);
    return {error, dError};
}

template<typename Scalar, typename UnderEstimation>
Scalar binarySearch(const PeriodicMatslise<Scalar> *ms, Scalar eMin, Scalar eMax, UnderEstimation under) {
    while (eMax - eMin > 1e-10) {
        Scalar e = (eMin + eMax) / 2;
        if (under(ms->matchingError(e))) {
            eMin = e;
        } else {
            eMax = e;
        }
    }

    return (eMin + eMax) / 2;
}

template<typename Scalar>
Scalar withinInterval(const PeriodicMatslise<Scalar> *ms, Scalar eMin, Scalar eMax) {
    Scalar fMin, fMax;
    fMin = calculateMatchingError(ms->matchingError(eMin).first).first;
    fMax = calculateMatchingError(ms->matchingError(eMax).first).first;
    if (fMin * fMax >= 0) {
        return abs(fMin) < abs(fMax) ? eMin : eMax;
    }

    bool segment = true;
    while (eMax - eMin > 1e-6) {
        Scalar e = segment ? (fMax * eMin - fMin * eMax) / (fMax - fMin) : (eMin + eMax) / 2;
        segment = !segment;
        Scalar f = calculateMatchingError(ms->matchingError(e).first).first;
        if (f == 0)
            return e;
        if (f * fMin < 0) {
            eMax = e;
            fMax = f;
        } else {
            eMin = e;
            fMin = f;
        }
    }

    return (eMin + eMax) / 2;
}

template<typename Scalar, bool with_eigenfunctions>
struct EigenpairsReturn {
    using type = vector<tuple<int, Scalar, std::conditional_t<with_eigenfunctions, vector<unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>, int>>>;
};

template<typename Scalar, bool with_eigenfunctions>
typename EigenpairsReturn<Scalar, with_eigenfunctions>::type
eigenpairsHelper(const PeriodicMatslise<Scalar> *ms, Scalar eMin, Scalar eMax, int iMin, int iMax) {
    Y<Scalar, 1, 2> boundary = Y<Scalar, 1, 2>::Periodic();
    typename EigenpairsReturn<Scalar, with_eigenfunctions>::type result;
    Scalar nextLow = binarySearch<Scalar>(ms, eMin, eMax, [iMin](const pair<Y<Scalar, 1, 2>, Array<Scalar, 2, 1>> &t) {
        return t.second.minCoeff() < iMin;
    });
    for (int i = iMin; i < iMax; ++i) {
        Scalar eLow = nextLow;
        Scalar eHigh = binarySearch<Scalar>(ms, eLow, eMax, [i](const pair<Y<Scalar, 1, 2>, Array<Scalar, 2, 1>> &t) {
            return t.second.maxCoeff() < i + 1;
        });
        // cout << "Searching between " << eLow << " and " << eHigh << endl;
        Scalar e = withinInterval<Scalar>(ms, eLow, eHigh);
        nextLow = binarySearch<Scalar>(ms, eHigh, eMax, [i](const pair<Y<Scalar, 1, 2>, Array<Scalar, 2, 1>> &t) {
            return t.second.minCoeff() < i + 1;
        });

        bool isDouble = i % 2 == 1 && nextLow - e < 1e-8;
        if constexpr (with_eigenfunctions) {
            auto &eigenfunctions = get<2>(
                    result.emplace_back(i, e, vector<unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>()));

            if (isDouble) {
                // double eigenvalue
                eigenfunctions.reserve(2);
                eigenfunctions.emplace_back(ms->matslise.eigenfunction(e, boundary.col(0), boundary.col(0)));
                eigenfunctions.emplace_back(ms->matslise.eigenfunction(e, boundary.col(1), boundary.col(1)));
            } else {
                EigenSolver<Matrix<Scalar, 2, 2>> solver(ms->matchingError(e).first.y(), true);
                Eigen::Index kernelIndex;
                solver.eigenvalues().array().abs().minCoeff(&kernelIndex);

                Y<Scalar> y = boundary * (Matrix<Scalar, 2, 1>) solver.eigenvectors().real().col(kernelIndex);

                eigenfunctions.reserve(1);
                eigenfunctions.emplace_back(ms->matslise.eigenfunction(e, y, y));
            }
        } else {
            result.emplace_back(i, e, isDouble ? 2 : 1);
        }
        if (isDouble) ++i;
    }
    return result;
}

template<typename Scalar, bool with_eigenfunctions>
typename EigenpairsReturn<Scalar, with_eigenfunctions>::type
eigenpairsByIndexHelper(const PeriodicMatslise<Scalar> *ms, int iMin, int iMax) {
    Scalar eMin = -1;
    Scalar eMax = 1;

    while (ms->matchingError(eMin).second.minCoeff() > iMin)
        eMin *= 2;

    while (ms->matchingError(eMax).second.maxCoeff() < iMax + 1)
        eMax *= 2;
    return eigenpairsHelper<Scalar, with_eigenfunctions>(ms, eMin, eMax, iMin, iMax);

}

template<typename Scalar>
std::vector<std::tuple<int, Scalar, std::vector<std::unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>>>
PeriodicMatslise<Scalar>::eigenpairsByIndex(int iMin, int iMax) const {
    return eigenpairsByIndexHelper<Scalar, true>(this, iMin, iMax);
}

template<typename Scalar>
std::vector<std::tuple<int, Scalar, int>> PeriodicMatslise<Scalar>::eigenvaluesByIndex(int iMin, int iMax) const {
    return eigenpairsByIndexHelper<Scalar, false>(this, iMin, iMax);
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
    cout << "eigenvalues: " << eigenvalues.transpose() << endl;
    cout << eigenvectors << endl;
    cout << endl;
    result.reserve(2);
    if (norm(eigenvalues[0]) < 1e-4) {
        Y<Scalar> boundary;
        boundary.y() = eigenvectors.col(0).real();
        result.emplace_back(matslise.eigenfunction(E, boundary, boundary));
    }
    if (norm(eigenvalues[1]) < 1e-4) {
        Y<Scalar> boundary;
        boundary.y() = eigenvectors.col(1).real();
        result.emplace_back(matslise.eigenfunction(E, boundary, boundary));
    }

    return result;
}

#include "instantiate.h"
