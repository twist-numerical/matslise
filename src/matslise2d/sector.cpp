#include <iostream>
#include <complex>
#include <queue>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/legendre.h"
#include "../util/horner.h"
#include "../util/constants.h"
#include "../util/calculateEta.h"
#include "./basisQuadrature.h"

using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace quadrature;


template<typename Scalar>
vector<typename Matscs<Scalar>::Sector> initializeMatscs(const typename Matslise2D<Scalar>::Sector &sector) {
    vector<typename Matscs<Scalar>::Sector> matscs;
    Index steps = sector.se2d->config.stepsPerSector;
    Scalar h = (sector.max - sector.min) / steps;
    matscs.reserve(steps);
    for (Index i = 0; i < steps; ++i) {
        Scalar min = sector.min + i * h;
        Scalar max = sector.max - (steps - i - 1) * h;

        matscs.emplace_back(
                legendre::getCoefficients<MATSCS_N, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Scalar>(
                        [&](Scalar y) -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> {
                            return sector.quadratures->dV(sector, y);
                        }, min, max),
                min, max, sector.direction);
    }
    return matscs;
}

template<typename Scalar>
Matslise2DSector<Scalar>::Matslise2DSector(const Matslise2D<Scalar> *se2d, const Scalar &ymin, const Scalar &ymax,
                                           Direction direction)
        : MatsliseNDSector<Scalar>(ymin, ymax), se2d(se2d), direction(direction) {
    MATSLISE_SCOPED_TIMER("2D sector");
    // std::cout << "new sector 2d" << std::endl;
    ybar = (ymax + ymin) / 2;
    function<Scalar(Scalar)> vbar_fun = [se2d, this](Scalar x) -> Scalar { return se2d->potential(x, ybar); };

    if (se2d->config.xSymmetric) {
        matslise = std::make_shared<MatsliseHalf<Scalar>>(
                vbar_fun, se2d->domain.template max<0>(), se2d->config.tolerance,
                sector_builder::getOrAutomatic<Matslise<Scalar>, false>(
                        se2d->config.xSectorBuilder, se2d->config.tolerance));
        quadratures = std::make_shared<BasisQuadrature<Scalar, MATSLISE2D_DELTA_V_DEGREE, true>>(
                static_cast<const MatsliseHalf<Scalar> *>(matslise.get())->ms);
    } else {
        matslise = std::make_shared<Matslise<Scalar>>(
                vbar_fun, se2d->domain.template min<0>(), se2d->domain.template max<0>(), se2d->config.tolerance,
                sector_builder::getOrAutomatic<Matslise<Scalar>, false>(
                        se2d->config.xSectorBuilder, se2d->config.tolerance));
        quadratures = std::make_shared<BasisQuadrature<Scalar, MATSLISE2D_DELTA_V_DEGREE, false>>(
                static_cast<const Matslise<Scalar> *>(matslise.get()));
    }

    const Index &N = se2d->config.basisSize;
    vector<pair<int, Scalar>> index_eigv = matslise->eigenvaluesByIndex(0, N, Y<Scalar>::Dirichlet());
    if (static_cast<int>(index_eigv.size()) != N) {
        throw std::runtime_error("SE2D: not enough basis-functions found on a sector");
    }
    eigenvalues.resize(N);
    eigenfunctions.resize(N);
    Scalar E;
    int index;
    for (int i = 0; i < N; ++i) {
        tie(index, E) = index_eigv[static_cast<unsigned long>(i)];
        eigenvalues[i] = E;
        // TODO: check i == index
        eigenfunctions[i] = matslise->eigenfunction(E, Y<Scalar>::Dirichlet(), index);
    }

    matscs = initializeMatscs<Scalar>(*this);
}

template<typename Scalar>
void Matslise2DSector<Scalar>::setDirection(Direction newDirection) {
    direction = newDirection;
    for (auto &sector : matscs)
        sector.setDirection(newDirection);
}

template<typename Scalar>
Matslise2DSector<Scalar> *Matslise2DSector<Scalar>::refine(
        const Matslise2D<Scalar> *problem, const Scalar &_min, const Scalar &_max, Direction _direction) const {
    Scalar h = _max - _min;
    if (ybar < _min + h / 3 || ybar > _max - h / 3) {
        return new Matslise2DSector<Scalar>(problem, _min, _max, _direction);
    }

    // std::cout << "refining sector 2d" << std::endl;
    auto sector = new Matslise2DSector<Scalar>(problem);
    sector->min = _min;
    sector->max = _max;
    sector->direction = _direction;
    sector->ybar = ybar;
    sector->matslise = matslise;
    sector->eigenfunctions = eigenfunctions;
    sector->eigenvalues = eigenvalues;
    sector->quadratures = quadratures;
    sector->matscs = initializeMatscs<Scalar>(*sector);
    return sector;
}

template<typename Scalar>
template<bool withDerivatives>
typename std::conditional<withDerivatives,
        std::pair<typename Matslise2DSector<Scalar>::ArrayXXs, typename Matslise2DSector<Scalar>::ArrayXXs>,
        typename Matslise2DSector<Scalar>::ArrayXXs>::type
Matslise2DSector<Scalar>::basis(const typename Matslise2DSector<Scalar>::ArrayXs &x) const {
    Eigen::Index size = x.size();

    const Index &N = se2d->config.basisSize;
    ArrayXXs b(size, N);
    ArrayXXs b_x;
    if constexpr(withDerivatives)
        b_x.resize(size, N);
    for (int i = 0; i < N; ++i) {
        Array<Y<Scalar>, Dynamic, 1> ys = eigenfunctions[i](x);
        b.col(i) = ys.template unaryExpr<std::function<Scalar(const Y<Scalar> &)>>(
                [](const Y<Scalar> &y) -> Scalar {
                    return y.data[0];
                });
        if constexpr(withDerivatives)
            b_x.col(i) = ys.template unaryExpr<std::function<Scalar(const Y<Scalar> &)>>(
                    [](const Y<Scalar> &y) -> Scalar {
                        return y.data[1];
                    });
    }
    if constexpr(withDerivatives)
        return {b, b_x};
    else
        return b;
}

template<typename Scalar>
template<bool withDerivatives>
typename std::conditional<withDerivatives,
        std::pair<typename Matslise2DSector<Scalar>::ArrayXs, typename Matslise2DSector<Scalar>::ArrayXs>,
        typename Matslise2DSector<Scalar>::ArrayXs>::type
Matslise2DSector<Scalar>::basis(const Scalar &x) const {
    ArrayXs b(this->eigenfunctions.size());
    ArrayXs b_x;
    if constexpr(withDerivatives)
        b_x.resize(this->eigenfunctions.size());

    const Index &N = se2d->config.basisSize;
    for (Index i = 0; i < N; ++i) {
        Y<Scalar> y = this->eigenfunctions[i](x);
        b[i] = y.data[0];
        if constexpr (withDerivatives)
            b_x[i] = y.data[1];
    }
    if constexpr (withDerivatives)
        return {b, b_x};
    else
        return b;
}

#define INSTANTIATE_BASIS(Scalar, r) \
template std::conditional<r, \
        std::pair<Matslise2DSector<Scalar>::ArrayXs, Matslise2DSector<Scalar>::ArrayXs>, \
        Matslise2DSector<Scalar>::ArrayXs>::type \
Matslise2DSector<Scalar>::basis<r>(const Scalar &) const; \
template std::conditional<r, \
        pair<Matslise2DSector<Scalar>::ArrayXXs, Matslise2DSector<Scalar>::ArrayXXs>, \
        Matslise2DSector<Scalar>::ArrayXXs \
    >::type \
Matslise2DSector<Scalar>::basis<r>(const Matslise2DSector<Scalar>::ArrayXs &) const;

#define INSTANTIATE_MORE(Scalar)\
INSTANTIATE_BASIS(Scalar, false)\
INSTANTIATE_BASIS(Scalar, true)

#include "../util/instantiate.h"