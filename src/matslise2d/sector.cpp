#include <iostream>
#include <complex>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/legendre.h"
#include "../util/calculateEta.h"
#include "../util/horner.h"
#include "./sector_correction_potential.h"

using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace quadrature;


template<typename Scalar>
Array<Scalar, MATSLISE_HMAX_delta, 2> getDeltaForEta(
        const typename Matslise<Scalar>::Sector *sector, const Scalar &E, const Y<Scalar> &y0) {
    Scalar Zd = sector->vs[0] - E;
    Array<Scalar, MATSLISE_HMAX_delta, MATSLISE_ETA_delta> u = sector->t_coeff.unaryExpr(
            [](const Matrix<Scalar, 2, 2> &m) -> Scalar { return m(0, 0); }).transpose();
    Array<Scalar, MATSLISE_HMAX_delta, MATSLISE_ETA_delta> v = sector->t_coeff.unaryExpr(
            [](const Matrix<Scalar, 2, 2> &m) -> Scalar { return m(0, 1); }).transpose();
    for (int i = MATSLISE_ETA_delta - 1; i >= 2; --i) {
        u.col(i - 2).template topRows<MATSLISE_HMAX_delta - 2>()
                += u.col(i).template bottomRows<MATSLISE_HMAX_delta - 2>() / Zd;
        u.col(i - 1).template topRows<MATSLISE_HMAX_delta - 2>()
                -= u.col(i).template bottomRows<MATSLISE_HMAX_delta - 2>() * (Scalar(2 * i - 3) / Zd);

        v.col(i - 2).template topRows<MATSLISE_HMAX_delta - 2>()
                += v.col(i).template bottomRows<MATSLISE_HMAX_delta - 2>() / Zd;
        v.col(i - 1).template topRows<MATSLISE_HMAX_delta - 2>()
                -= v.col(i).template bottomRows<MATSLISE_HMAX_delta - 2>() * (Scalar(2 * i - 3) / Zd);
    }



    Scalar h = sector->h*.33;
    Scalar *eta = calculateEta(Zd * h * h, 2);

    T<Scalar> t = sector->calculateT(E, h, true);
    cout << "\n --- " << endl;
    cout << "\n" << u.col(0).transpose() << endl;
    cout << u.col(1).transpose() << endl;
    cout << "\n" << v.col(0).transpose() << endl;
    cout << v.col(1).transpose() << endl;
    cout << "\n" << (
            horner<Scalar>(u.col(0), h, MATSLISE_HMAX_delta) * eta[0]
            + horner<Scalar>(u.col(1), h, MATSLISE_HMAX_delta) * eta[1]
    ) << ", " << (
                 horner<Scalar>(v.col(0), h, MATSLISE_HMAX_delta) * eta[0]
                 + horner<Scalar>(v.col(1), h, MATSLISE_HMAX_delta) * eta[1]
         ) << endl;
    cout << "\n" << t.t << endl;


    delete[] eta;

    return (y0.y(0) * u + y0.y(1) * v).template leftCols<2>();
};

template<typename Scalar>
Matslise2D<Scalar>::Sector::Sector(const Matslise2D<Scalar> *se2d, const Scalar &ymin, const Scalar &ymax,
                                   bool backward)
        : se2d(se2d), min(ymin), max(ymax), backward(backward) {
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
    vector<typename Matslise<Scalar>::Eigenfunction> func_eigenfunctions(se2d->N);
    Scalar E;
    int index;
    for (int i = 0; i < se2d->N; ++i) {
        tie(index, E) = index_eigv[static_cast<unsigned long>(i)];
        eigenvalues[i] = E;
        // TODO: check i == index
        func_eigenfunctions[i] = matslise->eigenfunction(E, Y<Scalar>::Dirichlet(), index);
        Array<Y<Scalar>, Dynamic, 1> func = func_eigenfunctions[i](se2d->grid);
        eigenfunctions[i] = ArrayXs(func.size());
        for (Eigen::Index j = 0; j < func.size(); ++j)
            eigenfunctions[i][j] = func[j].y[0];
    }

    if (auto m = dynamic_cast<Matslise<Scalar> *>(matslise)) {
        vector<vector<Array<Scalar, MATSLISE_HMAX_delta, 2>>> polynomials(m->sectors.size());
        auto polyIt = polynomials.begin();
        for (auto &sector : m->sectors) {
            polyIt->resize(se2d->N);
            auto nestedPolyIt = polyIt->begin();
            for (int i = 0; i < se2d->N; ++i) {
                if (!sector->backward) {
                    *nestedPolyIt = getDeltaForEta(sector, eigenvalues[i], func_eigenfunctions[i](sector->min));
                }
                ++nestedPolyIt;
            }
            ++polyIt;
        }

        polyIt = polynomials.begin();
        for (auto &sector : m->sectors) {
            if (!sector->backward) {
                ArrayXs grid = lobatto::grid<Scalar>(ArrayXs::LinSpaced(20, sector->min, sector->max));
                for (int i = 0; i < se2d->N; ++i)
                    for (int j = 0; j < se2d->N; ++j) {
                        Array<Scalar, MATSLISE_N, 1> quad = vbar_formulas((*polyIt)[i], (*polyIt)[j], sector->h,
                                                                          sector->vs[0] - eigenvalues[i],
                                                                          sector->vs[0] - eigenvalues[j]);

                        cout << "\n *** quad <-> lobatto" << endl;
                        cout << "\n" << quad.transpose() << endl;
                        Array<Scalar, Dynamic, 1> fi = func_eigenfunctions[i](grid).unaryExpr(
                                [](const Y<Scalar> &y) { return y.y[0]; });
                        Array<Scalar, Dynamic, 1> fj = func_eigenfunctions[j](grid).unaryExpr(
                                [](const Y<Scalar> &y) { return y.y[0]; });
                        cout << "\n" << lobatto::quadrature<Scalar>(grid, fi * fj) << endl;
                        cout << "\n" << lobatto::quadrature<Scalar>(grid, fi * fj * (grid - sector->min)) << endl;
                        cout << endl;
                    }

            }
            ++polyIt;
        }
    }

    matscs = new typename Matscs<Scalar>::Sector(
            legendre::getCoefficients<MATSCS_N, MatrixXs, Scalar>([this](Scalar y) -> MatrixXs {
                return this->calculateDeltaV(y);
            }, min, max),
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

    ArrayXs vDiff =
            se2d->grid.unaryExpr([this, y](const Scalar &x) -> Scalar { return this->se2d->potential(x, y); }) - vbar;

    for (int i = 0; i < se2d->N; ++i) {
        for (int j = 0; j <= i; ++j) {
            dV(i, j) = lobatto::quadrature<Scalar>(se2d->grid, eigenfunctions[i] * vDiff * eigenfunctions[j]);
            if (j < i) dV(j, i) = dV(i, j);
        }
        dV(i, i) += eigenvalues[i];
    }

    return dV;
/*
    std::function<ArrayXs(Scalar)> phi = this->template basis<false>();
    return trapezoidal::adaptive<Scalar, MatrixXs>([this, &phi, &y](const Scalar &x) -> MatrixXs {
                                                       ArrayXs f = phi(x);
                                                       return (se2d->potential(x, y) - matslise->potential(x)) * f.matrix() * f.matrix().transpose();
                                                   }, se2d->domain.sub.min, se2d->domain.sub.max, 1e-3,
                                                   [](const MatrixXs &m) -> Scalar { return m.cwiseAbs().maxCoeff(); });
                                                   */
}

template<typename Scalar>
template<bool withDerivative, typename diffType>
diffType Matslise2D<Scalar>::Sector::basis(const typename Matslise2D<Scalar>::ArrayXs &x) const {
    const Y<Scalar> y0 = Y<Scalar>({0, 1}, {0, 0});
    Eigen::Index size = x.size();

    ArrayXXs b(size, se2d->N);
    ArrayXXs b_x(size, se2d->N);
    for (int i = 0; i < se2d->N; ++i) {
        Array<Y<Scalar>, Dynamic, 1> ys = matslise->eigenfunction(eigenvalues[i], y0, i)(x);
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
        basis[static_cast<size_t>(index)] = matslise->eigenfunction(eigenvalues[index], y0, index);
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