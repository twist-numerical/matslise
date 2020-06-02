#include <iostream>
#include <complex>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/legendre.h"
#include "../util/horner.h"
#include "./sector_correction_potential.h"

using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace quadrature;


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
#define MATSLISE_N_integrate 12
#define MATSLISE_INTEGRATE_hmax 12
        vector<Array<Scalar, MATSLISE_N_integrate, 1>> quadData;
        quadData.reserve(se2d->N * (se2d->N + 1) / 2 * m->sectors.size());

        for (auto &sector : m->sectors) {
            Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> u = sector->t_coeff.unaryExpr(
                    [](const Matrix<Scalar, 2, 2, Eigen::DontAlign> &T) {
                        return T(0, 0);
                    });
            Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> v = sector->t_coeff.unaryExpr(
                    [](const Matrix<Scalar, 2, 2, Eigen::DontAlign> &T) {
                        return T(0, 1);
                    });

            Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_hmax, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> uu
                    = etaProduct<Scalar, MATSLISE_INTEGRATE_hmax>(u, u);
            Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_hmax, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> uv
                    = etaProduct<Scalar, MATSLISE_INTEGRATE_hmax>(u, v);
            Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_hmax, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> vv
                    = etaProduct<Scalar, MATSLISE_INTEGRATE_hmax>(v, v);

            std::vector<Y<Scalar>> y0;
            for (auto &f : func_eigenfunctions) {
                if (sector->backward) {
                    y0.emplace_back(f(sector->max));
                    y0.back().reverse();
                } else {
                    y0.emplace_back(f(sector->min));
                }
            }

           // ArrayXs grid = lobatto::grid<Scalar>(ArrayXs::LinSpaced(20, sector->min, sector->max));
            for (int i = 0; i < se2d->N; ++i)
                for (int j = 0; j <= i; ++j) {
                    Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta>
                            eta = (i == j
                                   ? eta_integrals<Scalar, true>(sector->h, sector->vs[0] - eigenvalues[i],
                                                                 sector->vs[0] - eigenvalues[j])
                                   : eta_integrals<Scalar, false>(sector->h, sector->vs[0] - eigenvalues[i],
                                                                  sector->vs[0] - eigenvalues[j]));

                    Eigen::Array<Scalar, MATSLISE_N_integrate, 1> quadratures = Eigen::Array<Scalar, MATSLISE_N_integrate, 1>::Zero();
                    {
                        Scalar suu = y0[i].y(0) * y0[j].y(0);
                        Scalar suv = y0[i].y(0) * y0[j].y(1);
                        Scalar svu = y0[i].y(1) * y0[j].y(0);
                        Scalar svv = y0[i].y(1) * y0[j].y(1);
                        for (int l = 0; l < MATSLISE_ETA_delta; ++l)
                            for (int m = 0; m < MATSLISE_ETA_delta; ++m) {
                                Eigen::Array<Scalar, MATSLISE_INTEGRATE_hmax, 1> s =
                                        (suu * uu(l, m) + suv * uv(l, m) + svu * uv(m, l) +
                                         svv * vv(l, m)).template head<MATSLISE_INTEGRATE_hmax>();
                                for (int n = std::max(0, 2 * l - 1 + 2 * m - 1); n < MATSLISE_INTEGRATE_hmax; ++n) {
                                    int count = std::min(MATSLISE_N_integrate, MATSLISE_INTEGRATE_delta - n);
                                    quadratures.segment(0, count) += s(n) * eta(l, m).segment(n, count);
                                }
                            }
                    }
                    legendreTransform<Scalar>(sector->h, quadratures);
                    quadData.push_back(quadratures);

                    if (sector->backward) {
                        auto &quad = quadData.back();
                        for (Eigen::Index k = 1; k < MATSLISE_N_integrate; k += 2)
                            quad(k) *= -1;
                    }

/*
                    auto &quad = quadData.back();
                    cout << "\n *** quad <-> lobatto" << endl;
                    cout << "\n" << quad.transpose() << endl;
                    Array<Scalar, Dynamic, 1> fi = func_eigenfunctions[i](grid).unaryExpr(
                            [](const Y<Scalar> &y) { return y.y[0]; });
                    Array<Scalar, Dynamic, 1> fj = func_eigenfunctions[j](grid).unaryExpr(
                            [](const Y<Scalar> &y) { return y.y[0]; });
                    cout << "\n" << lobatto::quadrature<Scalar>(grid, fi * fj) << endl;
                    cout << "\n"
                         << lobatto::quadrature<Scalar>(grid, fi * fj * (2 * (grid - sector->min) / sector->h - 1))
                         << endl;
                    cout << "\n"
                         << lobatto::quadrature<Scalar>(grid, fi * fj * (6 * (grid - sector->min) / sector->h *
                                                                         (grid - sector->min) / sector->h -
                                                                         6 * (grid - sector->min) / sector->h + 1))
                         << endl;

                    if (abs(lobatto::quadrature<Scalar>(grid, fi * fj) - quad(0)) > 1e-3)
                        cout << "wrong" << endl;
                    */

                }
        }

        matscs = new typename Matscs<Scalar>::Sector(
                legendre::getCoefficients<MATSCS_N, MatrixXs, Scalar>([&, this, m](Scalar y) -> MatrixXs {
                    MatrixXs dV = MatrixXs::Zero(this->se2d->N, this->se2d->N);
                    for (int i = 0; i < this->se2d->N; ++i) {
                        dV(i, i) = this->eigenvalues[i];
                    }

                    auto quadIt = quadData.begin();
                    for (auto &sector : m->sectors) {
                        //ArrayXs grid = lobatto::grid<Scalar>(ArrayXs::LinSpaced(30, sector->min, sector->max));
                        Array<Scalar, MATSLISE_N_integrate, 1> vDiff;

                        if (sector->backward) {
                            vDiff = Array<Scalar, MATSLISE_N_integrate, 1>::Map(sector->vs.data());
                            for (int i = 0; i < MATSLISE_N_integrate; i += 2)
                                vDiff(i) *= -1;
                        } else {
                            vDiff = -Array<Scalar, MATSLISE_N_integrate, 1>::Map(sector->vs.data());
                        }

                        vDiff += Array<Scalar, MATSLISE_N_integrate, 1>::Map(
                                legendre::getCoefficients<MATSLISE_N_integrate, Scalar, Scalar>(
                                        [this, y](const Scalar &x) -> Scalar { return this->se2d->potential(x, y); },
                                        sector->min, sector->max
                                ).data());

                        Scalar h = 1;
                        for (int i = 0; i < MATSLISE_N_integrate; ++i, h *= sector->h)
                            vDiff(i) *= h;

                        for (int i = 0; i < this->se2d->N; ++i)
                            for (int j = 0; j <= i; ++j) {
                                Scalar v = (quadIt->template head<MATSLISE_N_integrate>() * vDiff).sum();
                                dV(i, j) += v;
                                if (i != j)
                                    dV(j, i) += v;

/*
                                cout << "\n" << lobatto::quadrature<Scalar>(grid, grid.unaryExpr([&](const Scalar &x) {
                                    return func_eigenfunctions[i](x).y(0) * func_eigenfunctions[j](x).y(0) *
                                           (this->se2d->potential(x, y) - m->potential(x));
                                })) << endl;

                                cout << v << endl;

                                cout << sector->backward << endl;

                                cout << quadIt->transpose() << endl;

                                if (sector->backward)
                                    cout << "???" << endl;
                                */
                                ++quadIt;
                            }
                    }
/*
                    cout << "\n ---" << endl;
                    cout << "\n" << dV << endl;
                    cout << "\n" << this->calculateDeltaV(y) << endl;
*/                 // cout << "\n" << (dV - this->calculateDeltaV(y)) << endl << endl;

                    return dV;
                }, min, max),
                ymin, ymax, backward);
    } else {
        matscs = new typename Matscs<Scalar>::Sector(
                legendre::getCoefficients<MATSCS_N, MatrixXs, Scalar>([this](Scalar y) -> MatrixXs {
                    return this->calculateDeltaV(y);
                }, min, max),
                ymin, ymax, backward);
    }
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