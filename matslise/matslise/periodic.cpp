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

template<typename T>
bool same_sign(const T &a, const T &b) {
    return (a < 0 && b < 0) || (a > 0 && b > 0);
}

template<typename Scalar, typename Func, typename FuncRet = decltype(declval<Func>()(declval<Scalar>()))>
Scalar binarySearch(Scalar eMin, Scalar eMax, const Func &f, const Scalar &eps, FuncRet fMin) {
    Scalar mid = (eMin + eMax) / 2;
    while (eMax - eMin > eps) {
        FuncRet fMid = f(mid);
        if (same_sign(fMin, fMid)) {
            fMin = fMid;
            eMin = mid;
        } else
            eMax = mid;
        mid = (eMin + eMax) / 2;
    }

    return (eMin + eMax) / 2;
}

template<typename Scalar, typename Func>
Scalar binarySearch(Scalar eMin, Scalar eMax, const Func &f, const Scalar &eps) {
    return binarySearch<Scalar, Func>(eMin, eMax, f, eps, f(eMin));
}

template<typename Scalar>
pair<Scalar, Scalar> eigenvaluePair(const PeriodicMatslise<Scalar> *ms, Scalar Emin, Scalar Emax, Scalar eps) {
    Scalar h = 0.01 * (Emax - Emin);
    auto fError = [ms](Scalar E) {
        return calculateMatchingError(ms->matchingError(E).first).first;
    };
    Scalar minError = fError(Emin);
    Scalar fMin = (fError(Emin + h) - minError) / h;

    Scalar mid = binarySearch(Emin, Emax, [ms](Scalar E) {
        return calculateMatchingError(ms->matchingError(E).first).second;
    }, eps, fMin);

    Scalar midError = fError(mid);
    if (same_sign(minError, midError)) {
        if (abs(midError) > eps) cout << "PeriodicMatslise::eigenvaluePair wrong result" << endl;
        return {mid, mid};
    } else {
        return {binarySearch(Emin, mid, fError, eps, minError), binarySearch(mid, Emax, fError, eps, midError)};
    }
}

template<typename Scalar, bool with_eigenfunctions>
struct EigenpairsReturn {
    using type = vector<tuple<int, Scalar, std::conditional_t<with_eigenfunctions, vector<unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>, int>>>;
};

template<typename Scalar, bool with_eigenfunctions>
typename EigenpairsReturn<Scalar, with_eigenfunctions>::type
eigenpairsHelper(const PeriodicMatslise<Scalar> *ms, Scalar eMin, Scalar eMax, int iMin, int iMax, Scalar eps) {
    const Y<Scalar, 1, 2> boundary = Y<Scalar, 1, 2>::Periodic();
    typename EigenpairsReturn<Scalar, with_eigenfunctions>::type result;
    auto addToResult = [&](int i, Scalar E, bool isDouble) {
        if (i < iMin || i >= iMax || E < eMin || E >= eMax) return;
        if constexpr (with_eigenfunctions) {
            auto &eigenfunctions = get<2>(
                    result.emplace_back(i, E, vector<unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>()));

            if (isDouble) {
                // double eigenvalue
                eigenfunctions.reserve(2);
                eigenfunctions.emplace_back(ms->matslise.eigenfunction(E, boundary.col(0), boundary.col(0)));
                eigenfunctions.emplace_back(ms->matslise.eigenfunction(E, boundary.col(1), boundary.col(1)));
            } else {
                EigenSolver<Matrix<Scalar, 2, 2>> solver(ms->matchingError(E).first.y(), true);
                Eigen::Index kernelIndex;
                solver.eigenvalues().array().abs().minCoeff(&kernelIndex);

                Y<Scalar> y = boundary * (Matrix<Scalar, 2, 1>) solver.eigenvectors().real().col(kernelIndex);

                eigenfunctions.reserve(1);
                eigenfunctions.emplace_back(ms->matslise.eigenfunction(E, y, y));
            }
        } else {
            result.emplace_back(i, E, isDouble ? 2 : 1);
        }
    };

    int i = iMin - (iMin % 2);
    auto findLow = [&](Scalar E) {
        return ms->matchingError(E).second.minCoeff() - i + 1;
    };
    auto findHigh = [&](Scalar E) {
        return ms->matchingError(E).second.maxCoeff() - i - 1;
    };

    if (i == 0) {
        Scalar high = binarySearch<Scalar>(eMin, eMax, findHigh, eps);
        Scalar E = binarySearch(eMin, high, [ms](Scalar E) {
            return calculateMatchingError(ms->matchingError(E).first).first;
        }, eps);
        addToResult(0, E, false);
        i += 2;
    }
    Scalar low = eMin;
    for (; i <= iMax + 1; i += 2) {
        low = binarySearch<Scalar>(low, eMax, findLow, eps);
        Scalar high = binarySearch<Scalar>(low, eMax, findHigh, eps);
        pair<Scalar, Scalar> p = eigenvaluePair(ms, low, high, eps);
        if (p.second - p.first < eps) {
            addToResult(i - 1, (p.first + p.second) / 2, true);
        } else {
            addToResult(i - 1, p.first, false);
            addToResult(i, p.second, false);
        }
    }
    return result;
}

template<typename Scalar, bool with_eigenfunctions>
typename EigenpairsReturn<Scalar, with_eigenfunctions>::type
eigenpairsByIndexHelper(const PeriodicMatslise<Scalar> *ms, int iMin, int iMax) {
    Scalar eMin = -1;
    Scalar eMax = 1;

    while (ms->matchingError(eMin).second.minCoeff() > std::max(iMin - 1, 0))
        eMin *= 2;

    while (ms->matchingError(eMax).second.maxCoeff() < iMax + 2)
        eMax *= 2;

    return eigenpairsHelper<Scalar, with_eigenfunctions>(ms, eMin, eMax, iMin, iMax, 1e-8);
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