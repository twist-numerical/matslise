#include <iostream>
#include "../matslise.h"
#include "../util/lobatto.h"
#include "../util/legendre.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
Matslise2D<Scalar>::Sector::Sector(const Matslise2D<Scalar> *se2d, const Scalar &ymin, const Scalar &ymax,
                                   bool backward)
        : se2d(se2d), min(ymin), max(ymax) {
    const Scalar ybar = (ymax + ymin) / 2;
    function<Scalar(Scalar)> vbar_fun = [se2d, ybar](Scalar x) -> Scalar { return se2d->potential(x, ybar); };
    vbar = se2d->grid.unaryExpr(vbar_fun);
    if (se2d->options.nestedOptions._symmetric)
        matslise = new MatsliseHalf<Scalar>(vbar_fun, se2d->domain.sub.max, 1e-9, se2d->options.nestedOptions._builder);
    else
        matslise = new Matslise<Scalar>(vbar_fun, se2d->domain.sub.min, se2d->domain.sub.max, 1e-9,
                                        se2d->options.nestedOptions._builder);

    vector<pair<int, Scalar>> index_eigv = matslise->eigenvaluesByIndex(0, se2d->N, Y<Scalar>::Dirichlet());
    if (static_cast<int>(index_eigv.size()) != se2d->N) {
        throw std::runtime_error("SE2D: not enough basis-functions found on a sector");
    }
    eigenvalues = new Scalar[se2d->N];
    eigenfunctions = new ArrayXs[se2d->N];
    Scalar E;
    int index;
    for (int i = 0; i < se2d->N; ++i) {
        tie(index, E) = index_eigv[static_cast<unsigned long>(i)];
        eigenvalues[i] = E;
        // TODO: check i == index
        Array<Y<Scalar>, Dynamic, 1> func = matslise->eigenfunction(
                E, Y<Scalar>::Dirichlet(), se2d->grid, index);
        eigenfunctions[i] = ArrayXs(func.size());
        for (Eigen::Index j = 0; j < func.size(); ++j)
            eigenfunctions[i][j] = func[j].y[0];
    }

    matscs = new typename Matscs<Scalar>::Sector(
            legendre::getCoefficients<MATSCS_N, MatrixXs, Scalar>([this](Scalar y) -> MatrixXs {
                return this->calculateDeltaV(y); }, min, max),
                ymin, ymax, backward);
}

template<typename Scalar>
Matslise2D<Scalar>::Sector::~Sector() {
    delete matslise;
    delete matscs;
    delete[] eigenvalues;
    delete[] eigenfunctions;
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
typename Matslise2D<Scalar>::MatrixXs Matslise2D<Scalar>::Sector::calculateDeltaV(const Scalar &y) const {
    MatrixXs dV(se2d->N, se2d->N);

    ArrayXs vDiff = se2d->grid.unaryExpr([this, y](const Scalar &x) -> Scalar { return this->se2d->potential(x, y); }) - vbar;

    for (int i = 0; i < se2d->N; ++i) {
        for (int j = 0; j <= i; ++j) {
            dV(i, j) = lobatto::quadrature<Scalar>(se2d->grid, eigenfunctions[i] * vDiff * eigenfunctions[j]);
            if (j < i) dV(j, i) = dV(i, j);
        }
        dV(i, i) += eigenvalues[i];
    }

    return dV;
}

template<typename Scalar>
template<bool withDerivative, typename diffType>
diffType Matslise2D<Scalar>::Sector::basis(const typename Matslise2D<Scalar>::ArrayXs &x) const {
    const Y<Scalar> y0 = Y<Scalar>({0, 1}, {0, 0});
    Eigen::Index size = x.size();

    ArrayXXs b(size, se2d->N);
    ArrayXXs b_x(size, se2d->N);
    for (int i = 0; i < se2d->N; ++i) {
        auto ys = matslise->eigenfunction(eigenvalues[i], y0, x, i);
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
    const Y<Scalar> y0 = Y<Scalar>::Dirichlet(1);
    vector<function<Y<Scalar>(Scalar)>> basis(static_cast<size_t>(se2d->N));
    for (int index = 0; index < se2d->N; ++index) {
        basis[static_cast<size_t>(index)] = matslise->eigenfunctionCalculator(eigenvalues[index], y0, index);
    }
    return [basis](const Scalar &x) -> diffType {
        ArrayXs b(basis.size());
        ArrayXs b_x(basis.size());

        for (int i = 0; i < static_cast<int>(basis.size()); ++i) {
            Y<Scalar> y = basis[static_cast<size_t>(i)](x);
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