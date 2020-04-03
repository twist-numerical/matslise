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
    neumannBoundary = Y<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Neumann(se2d->N);
    dirichletBoundary = Y<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Dirichlet(se2d->N);
}

template<typename Scalar>
Matslise2DHalf<Scalar>::~Matslise2DHalf() {
    delete se2d;
}

template<typename Scalar>
Scalar Matslise2DHalf<Scalar>::firstEigenvalue() const {
    return se2d->firstEigenvalue(neumannBoundary);
}

template<typename Scalar>
Scalar Matslise2DHalf<Scalar>::eigenvalue(const Scalar &guessE) const {
    Scalar evenE = se2d->eigenvalue(neumannBoundary, guessE);
    Scalar oddE = se2d->eigenvalue(dirichletBoundary, guessE);
    if (abs(evenE - guessE) < abs(oddE - guessE))
        return evenE;
    else
        return oddE;
}

template<typename Scalar>
vector<Scalar> Matslise2DHalf<Scalar>::eigenvalues(const Scalar &Emin, const Scalar &Emax) const {
    vector<Scalar> eigenvalues = se2d->eigenvalues(neumannBoundary, Emin, Emax);
    vector<Scalar> eigenvaluesOdd = se2d->eigenvalues(dirichletBoundary, Emin, Emax);
    eigenvalues.insert(eigenvalues.end(),
                       make_move_iterator(eigenvaluesOdd.begin()),
                       make_move_iterator(eigenvaluesOdd.end()));
    sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
}

template<typename Scalar>
vector<typename Matslise2DHalf<Scalar>::ArrayXXs>
Matslise2DHalf<Scalar>::eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const {
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
        const Y<Scalar, Dynamic> &boundary = even ? neumannBoundary : dirichletBoundary;
        if (abs(E - se2d->eigenvalue(boundary, E)) < 1e-3) {
            vector<ArrayXXs> sortedFs = se2d->eigenfunction(boundary, E, x, sortedY);
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
vector<Scalar> Matslise2DHalf<Scalar>::eigenvaluesByIndex(int Imin, int Imax) const {
    vector<Scalar> valuesNeumann = se2d->eigenvaluesByIndex(neumannBoundary, 0, Imax);
    vector<Scalar> valuesDirichlet = se2d->eigenvaluesByIndex(dirichletBoundary, 0, Imax);
    vector<Scalar> result;
    merge(valuesNeumann.begin(), valuesNeumann.end(),
          valuesDirichlet.begin(), valuesDirichlet.end(),
          back_inserter(result));
    result.erase(result.begin() + Imax, result.end());
    result.erase(result.begin(), result.begin() + Imin);
    return result;
}

template<typename Scalar>
vector<function<Scalar(Scalar, Scalar)>> Matslise2DHalf<Scalar>::eigenfunctionCalculator(const Scalar &E) const {
    const Scalar SQRT1_2 = sqrt(Scalar(.5));
    vector<function<Scalar(Scalar, Scalar)>> result;
    for (bool even : {false, true}) {
        const Y<Scalar, Dynamic> &boundary = even ? neumannBoundary : dirichletBoundary;
        if (abs(E - se2d->eigenvalue(boundary, E)) < 1e-5) {
            for (const function<Scalar(Scalar, Scalar)> &f : se2d->eigenfunctionCalculator(boundary, E)) {
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