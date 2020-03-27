#include <iostream>
#include <map>
#include "../matslise.h"
#include "../util/lobatto.h"

using namespace Eigen;
using namespace matslise;
using namespace std;


template<typename Scalar>
Matslise2DHalf<Scalar>::Matslise2DHalf(const function<Scalar(const Scalar &, const Scalar &)> &V,
                                       const matslise::Rectangle<2, Scalar> &domain,
                                       const Options2<Scalar> &options) : domain(domain) {
    se2d = new Matslise2D<Scalar>(V, {domain.sub, 0, domain.getMax(1)}, options);
    evenBoundary = Y<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Neumann(se2d->N);
    oddBoundary = Y<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Dirichlet(se2d->N);
}

template<typename Scalar>
Matslise2DHalf<Scalar>::~Matslise2DHalf() {
    delete se2d;
}

template<typename Scalar>
void Matslise2DHalf<Scalar>::setParity(bool even) {
    se2d->y0Left = even ? evenBoundary : oddBoundary;
}

template<typename Scalar>
pair<Scalar, Scalar> Matslise2DHalf<Scalar>::matchingError(const Scalar &E) {
    setParity(true);
    pair<Scalar, Scalar> errorEven = se2d->matchingError(E);
    setParity(false);
    pair<Scalar, Scalar> errorOdd = se2d->matchingError(E);
    return abs(get<0>(errorEven)) <= abs(get<0>(errorOdd)) ? errorEven : errorOdd;
}

template<typename Scalar>
Scalar Matslise2DHalf<Scalar>::eigenvalue(
        const Scalar &guessE, const Scalar &tolerance, int maxIterations, const Scalar &minTolerance) {
    setParity(true);
    Scalar evenE = se2d->eigenvalue(guessE, tolerance, maxIterations, minTolerance);
    setParity(false);
    Scalar oddE = se2d->eigenvalue(guessE, tolerance, maxIterations, minTolerance);
    if (abs(evenE - guessE) < abs(oddE - guessE))
        return evenE;
    else
        return oddE;
}

template<typename Scalar>
vector<Scalar> Matslise2DHalf<Scalar>::eigenvalues(
        const Scalar &Emin, const Scalar &Emax, const int &initialSteps) {
    setParity(true);
    vector<Scalar> eigenvalues = se2d->eigenvalues(Emin, Emax, initialSteps);
    setParity(false);
    vector<Scalar> eigenvaluesOdd = se2d->eigenvalues(Emin, Emax, initialSteps);
    eigenvalues.insert(eigenvalues.end(),
                       make_move_iterator(eigenvaluesOdd.begin()),
                       make_move_iterator(eigenvaluesOdd.end()));
    sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
}

template<typename Scalar>
Scalar Matslise2DHalf<Scalar>::firstEigenvalue() {
    setParity(true);
    return se2d->firstEigenvalue();
}

template<typename Scalar>
vector<typename Matslise2DHalf<Scalar>::ArrayXXs>
Matslise2DHalf<Scalar>::eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) {
    Eigen::Index n = y.size();
    for (Eigen::Index i = 1; i < n; ++i)
        if (y[i - 1] > y[i])
            throw runtime_error("SE2DHalf::computeEigenfunction(): y has to be sorted");

    Array<bool, Dynamic, 1> isNegative(n);
    ArrayXs sortedY(n);
    {
        long i = n;
        long negativeIndex = 0;
        long positiveIndex = n;
        while (negativeIndex < positiveIndex) {
            --i;
            if (-y[negativeIndex] > y[positiveIndex - 1]) {
                isNegative[i] = true;
                sortedY[i] = -y[negativeIndex];
                negativeIndex += 1;
            } else {
                isNegative[i] = false;
                positiveIndex -= 1;
                sortedY[i] = y[positiveIndex];
            }
        }
    }

    const Scalar SQRT1_2 = sqrt(Scalar(.5));

    vector<ArrayXXs> result;
    for (bool even : {false, true}) {
        setParity(even);
        if (abs(E - se2d->eigenvalue(E)) < 1e-3) {
            vector<ArrayXXs> sortedFs = se2d->eigenfunction(E, x, sortedY);
            for (auto sortedF : sortedFs) {
                ArrayXXs f = ArrayXXs::Zero(x.size(), y.size());
                Eigen::Index negativeIndex = 0;
                Eigen::Index positiveIndex = n;
                for (Eigen::Index i = n - 1; i >= 0; --i) {
                    if (isNegative[i]) {
                        if (even)
                            f.col(negativeIndex) = sortedF.col(i);
                        else
                            f.col(negativeIndex) = -sortedF.col(i);
                        negativeIndex += 1;
                    } else {
                        positiveIndex -= 1;
                        f.col(positiveIndex) = sortedF.col(i);
                    }
                }
                f *= SQRT1_2;
                result.push_back(f);
            }

        }
    }
    return result;
}

template<typename Scalar>
vector<Scalar> Matslise2DHalf<Scalar>::eigenvaluesByIndex(int, int) {
    throw std::logic_error("Function not yet implemented");
}
// std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

template<typename Scalar>
vector<function<Scalar(Scalar, Scalar)>> Matslise2DHalf<Scalar>::eigenfunctionCalculator(const Scalar &E) {
    const Scalar SQRT1_2 = sqrt(Scalar(.5));
    vector<function<Scalar(Scalar, Scalar)>> result;
    for (bool even : {false, true}) {
        setParity(even);
        if (abs(E - se2d->eigenvalue(E)) < 1e-5) {
            for (const function<Scalar(Scalar, Scalar)> &f : se2d->eigenfunctionCalculator(E)) {
                result.push_back([f, even, SQRT1_2](const Scalar &x, const Scalar &y) -> Scalar {
                    Scalar r = f(x, y < 0 ? -y : y);
                    if (y < 0 && !even)
                        r = -r;
                    r *= SQRT1_2;
                    return r;
                });
            }
        }
    }
    return result;
}


#include "../util/instantiate.h"