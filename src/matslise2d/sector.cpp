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
typename Matscs<Scalar>::Sector *initializeMatscs(const typename Matslise2D<Scalar>::Sector &sector) {
    return new typename Matscs<Scalar>::Sector(
            legendre::getCoefficients<MATSCS_N, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Scalar>(
                    [&](Scalar y) -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> {
                        return sector.quadratures->dV(sector, y);
                    }, sector.min, sector.max),
            sector.min, sector.max, sector.direction);
}

template<typename Scalar>
Matslise2D<Scalar>::Sector::Sector(const Matslise2D<Scalar> *se2d, const Scalar &ymin, const Scalar &ymax,
                                   Direction direction)
        : se2d(se2d), min(ymin), max(ymax), direction(direction) {
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
void Matslise2D<Scalar>::Sector::setDirection(Direction newDirection) {
    direction = newDirection;
    matscs->setDirection(newDirection);
}

template<typename Scalar>
typename Matslise2D<Scalar>::Sector *Matslise2D<Scalar>::Sector::refine(
        const Matslise2D<Scalar> *problem, const Scalar &_min, const Scalar &_max, Direction _direction) const {
    Scalar h = _max - _min;
    if (ybar < _min + h / 3 || ybar > _max - h / 3) {
        return new Sector(problem, _min, _max, _direction);
    }

    // std::cout << "refining sector 2d" << std::endl;
    auto sector = new Sector(problem);
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
void clamp(Scalar &value, const Scalar &min, const Scalar &max) {
    if (value < min)
        value = min;
    else if (value > max)
        value = max;
}


template<typename Scalar, typename Derived>
Matrix<complex<Scalar>, Dynamic, Dynamic> theta(
        const MatrixBase<Derived> &U,
        const MatrixBase<Derived> &V) {
    return (V - U * complex<Scalar>(0, 1))
            .transpose()
            .partialPivLu()
            .solve((V + U * complex<Scalar>(0, 1)).transpose())
            .transpose();
}

template<typename Scalar>
inline Array<Scalar, Dynamic, 1> angle(const Matrix<complex<Scalar>, Dynamic, Dynamic> &m) {
    const Scalar PI2 = constants<Scalar>::PI * 2;
    return m.eigenvalues().array().arg().unaryExpr(
            [&](const Scalar &a) -> Scalar { return a < -1e-16 ? a + PI2 : a; });
}

template<typename Scalar>
Index Matslise2D<Scalar>::Sector::estimateIndex(
        const Scalar &E, const Y<Scalar, Eigen::Dynamic> &y0, const Y<Scalar, Eigen::Dynamic> &y1) const {
    Index n = matscs->n;
    using ArrayXs = typename Matslise2D<Scalar>::ArrayXs;
    using MatrixXcs = Matrix<complex<Scalar>, Dynamic, Dynamic>;

    Index zeros = 0;
    queue<tuple<int, Scalar, Scalar, MatrixXs, MatrixXs, ArrayXs>> todo;
    todo.emplace(0, min, max, y0.getY(0), y0.getY(1), angle<Scalar>(theta<Scalar>(y1.getY(0), y1.getY(1))));

    while (!todo.empty()) {
        const int &depth = get<0>(todo.front());
        const Scalar &a = get<1>(todo.front());
        const Scalar &b = get<2>(todo.front());
        Scalar h = b - a;
        const MatrixXs &U0 = get<3>(todo.front());
        const MatrixXs &V0 = get<4>(todo.front());
        const ArrayXs &betas = get<5>(todo.front());

        ArrayXs Z = h * h * (matscs->vs[0].diagonal().array() - ArrayXs::Constant(n, E));
        Array<Scalar, 2, Dynamic> eta = calculateEta<Scalar, 2>(Z);

        MatrixXcs thetaZ0 = theta<Scalar>(U0, V0);
        const complex<Scalar> i_delta(0, h);
        ArrayXs alphas = angle<Scalar>(
                ((eta.row(0) + i_delta * eta.row(1)) / (eta.row(0) - i_delta * eta.row(1))).matrix().asDiagonal() *
                thetaZ0
        );
        if (depth < 2 && ((betas < 1e-4).any() || ((betas - alphas).abs() > 6.283).any())) {
            // Zeroth order propagation probably inaccurate: refine steps
            Scalar mid = (a + b) / 2;
            Y<Scalar, Dynamic> yMid = direction == forward
                                      ? propagate(E, y0, a, mid)
                                      : propagate(E, y1, b, mid);
            todo.emplace(depth + 1, a, mid, U0, V0, angle<Scalar>(theta<Scalar>(yMid.getY(0), yMid.getY(1))));
            todo.emplace(depth + 1, mid, b, yMid.getY(0), yMid.getY(1), betas);
        } else {
            Scalar argdet = 0;
            for (int i = 0; i < n; ++i) {
                if (Z[i] < 0) {
                    Scalar sZ = sqrt(-Z[i]);
                    argdet += (sZ + atan2(
                            (h - sZ) * eta(1, i) * eta(0, i),
                            1 + (h * sZ + Z[i]) * eta(1, i) * eta(1, i)));
                } else {
                    argdet += atan2(h * eta(1, i), eta(0, i));
                }
            }
            argdet *= 2;

            zeros += round((angle<Scalar>(thetaZ0).sum() + argdet - alphas.sum()) / (2 * constants<Scalar>::PI));
        }
        todo.pop();
    }
    return zeros;
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
template<bool withDerivatives>
typename std::conditional<withDerivatives,
        std::pair<typename Matslise2D<Scalar>::ArrayXXs, typename Matslise2D<Scalar>::ArrayXXs>,
        typename Matslise2D<Scalar>::ArrayXXs>::type
Matslise2D<Scalar>::Sector::basis(const typename Matslise2D<Scalar>::ArrayXs &x) const {
    Eigen::Index size = x.size();

    ArrayXXs b(size, se2d->N);
    ArrayXXs b_x;
    if constexpr(withDerivatives)
        b_x.resize(size, se2d->N);
    for (int i = 0; i < se2d->N; ++i) {
        Array<Y<Scalar>, Dynamic, 1> ys = eigenfunctions[i](x);
        b.col(i) = ys.template unaryExpr<std::function<Scalar(const Y<Scalar> &)>>(
                [](const Y<Scalar> &y) -> Scalar {
                    return y.y[0];
                });
        if constexpr(withDerivatives)
            b_x.col(i) = ys.template unaryExpr<std::function<Scalar(const Y<Scalar> &)>>(
                    [](const Y<Scalar> &y) -> Scalar {
                        return y.y[1];
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
        std::pair<typename Matslise2D<Scalar>::ArrayXs, typename Matslise2D<Scalar>::ArrayXs>,
        typename Matslise2D<Scalar>::ArrayXs>::type
Matslise2D<Scalar>::Sector::basis(const Scalar &x) const {
    ArrayXs b(this->eigenfunctions.size());
    ArrayXs b_x;
    if constexpr(withDerivatives)
        b_x.resize(this->eigenfunctions.size());

    for (Index i = 0; i < this->se2d->N; ++i) {
        Y<Scalar> y = this->eigenfunctions[i](x);
        b[i] = y.y[0];
        if constexpr (withDerivatives)
            b_x[i] = y.y[1];
    }
    if constexpr (withDerivatives)
        return {b, b_x};
    else
        return b;
}

#define INSTANTIATE_PROPAGATE(Scalar, r) \
template Y<Scalar, Dynamic, r> \
Matslise2D<Scalar>::Sector::propagate<r>( \
        const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &, const Scalar &, const Scalar &, bool) const;

#define INSTANTIATE_BASIS(Scalar, r) \
template std::conditional<r, \
        std::pair<Matslise2D<Scalar>::ArrayXs, Matslise2D<Scalar>::ArrayXs>, \
        Matslise2D<Scalar>::ArrayXs>::type \
Matslise2D<Scalar>::Sector::basis<r>(const Scalar &) const; \
template std::conditional<r, \
        pair<Matslise2D<Scalar>::ArrayXXs, Matslise2D<Scalar>::ArrayXXs>, \
        Matslise2D<Scalar>::ArrayXXs \
    >::type \
Matslise2D<Scalar>::Sector::basis<r>(const Matslise2D<Scalar>::ArrayXs &) const;

#define INSTANTIATE_MORE(Scalar)\
INSTANTIATE_PROPAGATE(Scalar, 1)\
INSTANTIATE_PROPAGATE(Scalar, -1)\
INSTANTIATE_BASIS(Scalar, false)\
INSTANTIATE_BASIS(Scalar, true)

#include "../util/instantiate.h"