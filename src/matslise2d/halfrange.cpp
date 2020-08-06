#include <iostream>
#include <map>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/constants.h"

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
pair<Scalar, Index> Matslise2DHalf<Scalar>::eigenvalue(const Scalar &guessE) const {
    pair<Scalar, Index> evenE = se2d->eigenvalue(neumannBoundary, guessE);
    pair<Scalar, Index> oddE = se2d->eigenvalue(dirichletBoundary, guessE);
    if (abs(evenE.first - guessE) < abs(oddE.first - guessE))
        return evenE;
    else
        return oddE;
}

template<typename Scalar>
Scalar Matslise2DHalf<Scalar>::eigenvalueError(const Scalar &E) const {
    return std::min(se2d->eigenvalueError(neumannBoundary, E), se2d->eigenvalueError(dirichletBoundary, E));
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> Matslise2DHalf<Scalar>::eigenvalues(const Scalar &Emin, const Scalar &Emax) const {
    vector<tuple<Index, Scalar, Index>> eigenvalues = se2d->eigenvalues(neumannBoundary, Emin, Emax);
    vector<tuple<Index, Scalar, Index>> eigenvaluesOdd = se2d->eigenvalues(dirichletBoundary, Emin, Emax);
    eigenvalues.insert(eigenvalues.end(),
                       make_move_iterator(eigenvaluesOdd.begin()),
                       make_move_iterator(eigenvaluesOdd.end()));
    sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> Matslise2DHalf<Scalar>::eigenvaluesByIndex(Index Imin, Index Imax) const {
    vector<tuple<Index, Scalar, Index>> valuesNeumann = se2d->eigenvaluesByIndex(neumannBoundary, 0, Imax);
    vector<tuple<Index, Scalar, Index>> valuesDirichlet = se2d->eigenvaluesByIndex(dirichletBoundary, 0, Imax);
    vector<tuple<Index, Scalar, Index>> result;
    merge(valuesNeumann.begin(), valuesNeumann.end(),
          valuesDirichlet.begin(), valuesDirichlet.end(),
          back_inserter(result));
    result.erase(result.begin() + Imax, result.end());
    result.erase(result.begin(), result.begin() + Imin);
    return result;
}


template<typename Scalar>
template<bool withDerivatives>
vector<Eigenfunction2D<Scalar, withDerivatives>> Matslise2DHalf<Scalar>::eigenfunction(const Scalar &E) const {
    vector<Eigenfunction2D<Scalar, withDerivatives>> eigenfunctions;
    for (bool even : {false, true}) {
        const Y<Scalar, Dynamic> &boundary = even ? neumannBoundary : dirichletBoundary;
        if (abs(E - se2d->eigenvalue(boundary, E).first) < 1e-3) {
            for (const auto &f : se2d->template eigenfunction<withDerivatives>(boundary, E)) {
                eigenfunctions.push_back(
                        {
                                [even, f](const Scalar &x, const Scalar &y) {
                                    auto r = f(x, y < 0 ? -y : y);
                                    Scalar scale =
                                            y < 0 && !even ? -constants<Scalar>::SQRT1_2 : constants<Scalar>::SQRT1_2;
                                    if constexpr(withDerivatives) {
                                        get<0>(r) *= scale;
                                        get<1>(r) *= scale;
                                        get<2>(r) *= scale;
                                    } else {
                                        r *= scale;
                                    }
                                    return r;
                                },
                                [even, f](const ArrayXs &x, const ArrayXs &y)
                                        -> typename Eigenfunction2D<Scalar, withDerivatives>::ArrayReturn {
                                    typename Eigenfunction2D<Scalar, withDerivatives>::ArrayReturn result
                                            = f(x, y.abs());
                                    if constexpr(withDerivatives) {
                                        get<0>(result) *= constants<Scalar>::SQRT1_2;
                                        get<1>(result) *= constants<Scalar>::SQRT1_2;
                                        get<2>(result) *= constants<Scalar>::SQRT1_2;
                                    } else {
                                        result *= constants<Scalar>::SQRT1_2;
                                    }
                                    if (!even) {
                                        for (int i = 0; i < y.size(); ++i) {
                                            if (y[i] < 0) {
                                                if constexpr(withDerivatives) {
                                                    get<0>(result).col(i) *= -1;
                                                    get<1>(result).col(i) *= -1;
                                                    get<2>(result).col(i) *= -1;
                                                } else {
                                                    result.col(i) *= -1;
                                                }
                                            }
                                        }
                                    }
                                    return result;
                                }
                        });
            }
        }
    }
    return eigenfunctions;
}


#include "../util/instantiate.h"