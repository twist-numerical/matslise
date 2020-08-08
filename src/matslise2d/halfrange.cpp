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
Index Matslise2DHalf<Scalar>::estimateIndex(const Scalar &E) const {
    return se2d->estimateIndex(dirichletBoundary, E) + se2d->estimateIndex(neumannBoundary, E);
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> mergeEigenvalues(
        const Matslise2DHalf<Scalar> &problem, const Scalar &Emin,
        const vector<tuple<Index, Scalar, Index>> &neumann, const vector<tuple<Index, Scalar, Index>> &dirichlet) {
    vector<tuple<Index, Scalar, Index>> result;
    if (neumann.empty() && dirichlet.empty())
        return result;

    Index first;
    if (neumann.empty())
        first = get<0>(dirichlet.front()) + problem.se2d->estimateIndex(problem.neumannBoundary, Emin);
    else if (dirichlet.empty())
        first = get<0>(neumann.front()) + problem.se2d->estimateIndex(problem.dirichletBoundary, Emin);
    else
        first = get<0>(neumann.front()) + get<0>(dirichlet.front());

    auto even = neumann.begin();
    auto odd = dirichlet.begin();
    Index currentIndex = first;
    while (even != neumann.end() || odd != dirichlet.end()) {
        auto &smallest =
                odd == dirichlet.end() || (even != neumann.end() && get<1>(*even) < get<1>(*odd))
                ? even : odd;

        const Scalar &E = get<1>(*smallest);
        const Index &multiplicity = get<2>(*smallest);
        if (!result.empty() && abs(get<1>(result.back()) - E) < 1e-4) {
            get<2>(result.back()) += multiplicity;
        } else {
            result.emplace_back(currentIndex, E, multiplicity);
        }
        currentIndex += multiplicity;
        ++smallest;
    }

    return result;
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> Matslise2DHalf<Scalar>::eigenvalues(const Scalar &Emin, const Scalar &Emax) const {
    return mergeEigenvalues<Scalar>(
            *this, Emin,
            se2d->eigenvalues(neumannBoundary, Emin, Emax), se2d->eigenvalues(dirichletBoundary, Emin, Emax));
}

template<typename Scalar, bool upper>
Scalar findEigenvalueForIndex(const Matslise2DHalf<Scalar> &problem, Index i) {
    Scalar Emin = 0;
    Index Imin;
    Scalar step = 1;
    while ((Imin = problem.estimateIndex(Emin)) > i) {
        Emin -= step;
        step *= 2;
    }
    step = 1;
    Scalar Emax = Emin + step;
    Index Imax;
    while ((Imax = problem.estimateIndex(Emax)) <= i) {
        step *= 2;
        Emax += step;
    }
    while (Emax - Emin > 1e-4) {
        Scalar Emid = (Emax + Emin) / 2;
        Index Imid = problem.estimateIndex(Emid);
        if (Imid <= i) {
            Emin = Emid;
            Imin = Imid;
        } else {
            Emax = Emid;
            Imax = Imid;
        }
    }
    return upper ? Emax : Emin;
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> Matslise2DHalf<Scalar>::eigenvaluesByIndex(Index Imin, Index Imax) const {
    return eigenvalues(findEigenvalueForIndex<Scalar, false>(*this, Imin),
                       findEigenvalueForIndex<Scalar, true>(*this, Imax));
}

template<typename Scalar>
template<bool withDerivatives>
vector<Eigenfunction2D<Scalar, withDerivatives>> Matslise2DHalf<Scalar>::eigenfunction(const Scalar &E) const {
    vector<Eigenfunction2D<Scalar, withDerivatives>> eigenfunctions;
    for (bool even : {false, true}) {
        const Y<Scalar, Dynamic> &boundary = even ? neumannBoundary : dirichletBoundary;
        if (abs(E - se2d->eigenvalue(boundary, E).first) < 1e-4) {
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