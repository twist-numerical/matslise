#include <iostream>
#include <complex>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/legendre.h"
#include "../util/horner.h"
#include "./basisQuadrature.h"

using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace quadrature;


template<typename Scalar>
typename Matscs<Scalar>::Sector *initializeMatscs(const typename Matslise2D<Scalar>::Sector &sector) {
    return new typename Matscs<Scalar>::Sector(
            legendre::getCoefficients<MATSCS_N, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Scalar>(
                    [&](Scalar y) -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> {
                        return sector.quadratures->dV(sector, y);
                    }, sector.min, sector.max),
            sector.min, sector.max, sector.backward);
}

template<typename Scalar>
Matslise2D<Scalar>::Sector::Sector(const Matslise2D<Scalar> *se2d, const Scalar &ymin, const Scalar &ymax,
                                   bool backward)
        : se2d(se2d), min(ymin), max(ymax), backward(backward) {
    // std::cout << "new sector 2d" << std::endl;
    ybar = (ymax + ymin) / 2;
    function<Scalar(Scalar)> vbar_fun = [se2d, this](Scalar x) -> Scalar { return se2d->potential(x, ybar); };

    if (se2d->options.nestedOptions._symmetric) {
        matslise = std::make_shared<MatsliseHalf<Scalar>>(
                vbar_fun, se2d->domain.sub.max, 1e-9, se2d->options.nestedOptions._builder);
        quadratures = std::make_shared<BasisQuadrature<Scalar, 8, true>>(
                static_cast<const MatsliseHalf<Scalar> *>(matslise.get())->ms);
    } else {
        matslise = std::make_shared<Matslise<Scalar>>(
                vbar_fun, se2d->domain.sub.min, se2d->domain.sub.max, 1e-9, se2d->options.nestedOptions._builder);
        quadratures = std::make_shared<BasisQuadrature<Scalar, 8, false>>(
                static_cast<const Matslise<Scalar> *>(matslise.get()));
    }

    vector<pair<int, Scalar>> index_eigv = matslise->eigenvaluesByIndex(0, se2d->N, Y<Scalar>::Dirichlet());
    if (static_cast<int>(index_eigv.size()) != se2d->N) {
        throw std::runtime_error("SE2D: not enough basis-functions found on a sector");
    }
    eigenvalues.resize(se2d->N);
    eigenfunctions.resize(se2d->N);
    Scalar E;
    int index;
    for (int i = 0; i < se2d->N; ++i) {
        tie(index, E) = index_eigv[static_cast<unsigned long>(i)];
        eigenvalues[i] = E;
        // TODO: check i == index
        eigenfunctions[i] = matslise->eigenfunction(E, Y<Scalar>::Dirichlet(), index);
    }

    matscs = initializeMatscs<Scalar>(*this);
}

template<typename Scalar>
Matslise2D<Scalar>::Sector::~Sector() {
    delete matscs;
}

template<typename Scalar>
typename Matslise2D<Scalar>::Sector *Matslise2D<Scalar>::Sector::refine(
        const Matslise2D<Scalar> *problem, const Scalar &_min, const Scalar &_max, bool _backward) const {
    Scalar h = _max - _min;
    if (backward != _backward || ybar < _min + h / 3 || ybar > _max - h / 3) {
        return new Sector(problem, _min, _max, _backward);
    }

    // std::cout << "refining sector 2d" << std::endl;
    auto sector = new Sector(problem);
    sector->min = _min;
    sector->max = _max;
    sector->backward = _backward;
    sector->ybar = ybar;
    sector->matslise = matslise;
    sector->eigenfunctions = eigenfunctions;
    sector->eigenvalues = eigenvalues;
    sector->quadratures = quadratures;
    sector->matscs = initializeMatscs<Scalar>(*sector);
    return sector;
}

template<typename Scalar>
template<int r>
Y<Scalar, Eigen::Dynamic, r>
Matslise2D<Scalar>::Sector::propagate(
        const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b, bool use_h) const {
    return matscs->propagateColumn(
            E, y0, a < min ? min : a > max ? max : a, b < min ? min : b > max ? max : b, use_h);
}

template<typename Scalar>
Scalar Matslise2D<Scalar>::Sector::error() const {
    return matscs->error();
}

template<typename Scalar>
template<bool withDerivative, typename diffType>
diffType Matslise2D<Scalar>::Sector::basis(const typename Matslise2D<Scalar>::ArrayXs &x) const {
    Eigen::Index size = x.size();

    ArrayXXs b(size, se2d->N);
    ArrayXXs b_x(size, se2d->N);
    for (int i = 0; i < se2d->N; ++i) {
        Array<Y<Scalar>, Dynamic, 1> ys = eigenfunctions[i](x);
        b.col(i) = ys.template unaryExpr<std::function<Scalar(const Y<Scalar> &)>>(
                [](const Y<Scalar> &y) -> Scalar {
                    return y.y[0];
                });
        if constexpr(withDerivative)
            b_x.col(i) = ys.template unaryExpr<std::function<Scalar(const Y<Scalar> &)>>(
                    [](const Y<Scalar> &y) -> Scalar {
                        return y.y[1];
                    });
    }
    if constexpr(withDerivative)
        return {b, b_x};
    else
        return b;
}

template<typename Scalar>
template<bool withDerivative, typename diffType>
function<diffType(Scalar)> Matslise2D<Scalar>::Sector::basis() const {
    return [this](const Scalar &x) -> diffType {
        ArrayXs b(this->eigenfunctions.size());
        ArrayXs b_x(this->eigenfunctions.size());

        for (int i = 0; i < static_cast<int>(this->eigenfunctions.size()); ++i) {
            Y<Scalar> y = this->eigenfunctions[static_cast<size_t>(i)](x);
            b[i] = y.y[0];
            if constexpr (withDerivative)
                b_x[i] = y.y[1];
        }
        if constexpr (withDerivative)
            return {b, b_x};
        else
            return b;
    };
}

#define INSTANTIATE_PROPAGATE(Scalar, r) \
template Y<Scalar, Dynamic, r> \
Matslise2D<Scalar>::Sector::propagate<r>( \
        const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b, bool use_h) const;

#define INSTANTIATE_BASIS(Scalar, r) \
template function<std::conditional<r, \
        std::pair<Matslise2D<Scalar>::ArrayXs, Matslise2D<Scalar>::ArrayXs>, \
        Matslise2D<Scalar>::ArrayXs>::type(Scalar)> \
Matslise2D<Scalar>::Sector::basis<r, std::conditional<r, \
    std::pair<Matslise2D<Scalar>::ArrayXs, Matslise2D<Scalar>::ArrayXs>, \
    Matslise2D<Scalar>::ArrayXs>::type>() const; \
template std::conditional<r, \
        pair<Matslise2D<Scalar>::ArrayXXs, Matslise2D<Scalar>::ArrayXXs>, \
        Matslise2D<Scalar>::ArrayXXs \
    >::type \
Matslise2D<Scalar>::Sector::basis<r, std::conditional<r, \
        pair<Matslise2D<Scalar>::ArrayXXs, Matslise2D<Scalar>::ArrayXXs>, \
        Matslise2D<Scalar>::ArrayXXs \
    >::type>(const Matslise2D<Scalar>::ArrayXs &x) const;

#define INSTANTIATE_MORE(Scalar)\
INSTANTIATE_PROPAGATE(Scalar, 1)\
INSTANTIATE_PROPAGATE(Scalar, -1)\
INSTANTIATE_BASIS(Scalar, false)\
INSTANTIATE_BASIS(Scalar, true)

#include "../util/instantiate.h"