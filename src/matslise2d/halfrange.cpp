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
                                       const Options2<Scalar> &options) : AbstractMatslise2D<Scalar>(V, domain) {
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
Scalar Matslise2DHalf<Scalar>::eigenvalueError(const Scalar &E) const {
    return std::min(se2d->eigenvalueError(neumannBoundary, E), se2d->eigenvalueError(dirichletBoundary, E));
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
template<bool withDerivative, typename returnType>
returnType Matslise2DHalf<Scalar>::eigenfunctionHelper(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const {
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

    returnType result;
    for (bool even : {false, true}) {
        const Y<Scalar, Dynamic> &boundary = even ? neumannBoundary : dirichletBoundary;
        if (abs(E - se2d->eigenvalue(boundary, E)) < 1e-3) {
            auto sortedFs = se2d->template eigenfunctionHelper<withDerivative>(boundary, E, x, sortedY);
            for (auto sortedF : sortedFs) {
                typename std::conditional<withDerivative, std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>, ArrayXXs>::type f;
                if constexpr(withDerivative)
                    f = {ArrayXXs::Zero(x.size(), y.size()), ArrayXXs::Zero(x.size(), y.size()),
                         ArrayXXs::Zero(x.size(), y.size())};
                else
                    f = ArrayXXs::Zero(x.size(), y.size());
                Eigen::Index negativeIndex = 0;
                Eigen::Index positiveIndex = n;
                for (Eigen::Index i = n - 1; i >= 0; --i) {
                    Eigen::Index index = isNegative[i] ? negativeIndex++ : --positiveIndex;
                    Scalar scale = !isNegative[i] || even ? SQRT1_2 : -SQRT1_2;
                    if constexpr(withDerivative) {
                        get<0>(f).col(index) = scale * get<0>(sortedF).col(i);
                        get<1>(f).col(index) = scale * get<1>(sortedF).col(i);
                        get<2>(f).col(index) = scale * get<2>(sortedF).col(i);
                    } else {
                        f.col(index) = scale * sortedF.col(i);
                    }
                }
                result.push_back(f);
            }

        }
    }
    return result;
}

template<typename Scalar>
template<bool withDerivative, typename returnType>
returnType Matslise2DHalf<Scalar>::eigenfunctionHelper(const Scalar &E) const {
    const Scalar SQRT1_2 = sqrt(Scalar(.5));
    returnType result;
    for (bool even : {false, true}) {
        const Y<Scalar, Dynamic> &boundary = even ? neumannBoundary : dirichletBoundary;
        if (abs(E - se2d->eigenvalue(boundary, E)) < 1e-5) {
            for (const auto &f : se2d->template eigenfunctionHelper<withDerivative>(boundary, E)) {
                result.push_back([f, even, SQRT1_2](const Scalar &x, const Scalar &y) -> auto {
                    auto r = f(x, y < 0 ? -y : y);
                    Scalar scale = y < 0 && !even ? -SQRT1_2 : SQRT1_2;
                    if constexpr(withDerivative) {
                        get<0>(r) *= scale;
                        get<1>(r) *= scale;
                        get<2>(r) *= scale;
                    } else {
                        r *= scale;
                    }
                    return r;
                });
            }
        }
    }
    return result;
}


#include "../util/instantiate.h"