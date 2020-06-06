//
// Created by toon on 6/3/20.
//

#ifndef MATSLISE_BASISQUADRATURE_H
#define MATSLISE_BASISQUADRATURE_H


using namespace std;
using namespace Eigen;
using namespace matslise;
using namespace quadrature;

template<typename Scalar>
class AbstractBasisQuadrature {
public:
    virtual Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dV(const Scalar &y) = 0;
};

template<typename Scalar, int hmax, bool skipOdd = false>
class BasisQuadrature : public AbstractBasisQuadrature<Scalar> {
public:
    const typename matslise::Matslise2D<Scalar>::Sector *sector2d;
    const matslise::Matslise<Scalar> *matslise;
    std::vector<Eigen::Array<Scalar, hmax, 1>> quadData;

    BasisQuadrature() {}

    BasisQuadrature(
            const typename matslise::Matslise2D<Scalar>::Sector *sector2d,
            const matslise::Matslise<Scalar> *matslise) : sector2d(sector2d), matslise(matslise) {
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dV(const Scalar &y) override {
        if (quadData.empty())
            calculateQuadData();

        Eigen::Index N = sector2d->se2d->N;
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dV
                = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);

        for (int i = 0; i < N; ++i) {
            dV(i, i) = sector2d->eigenvalues[i];
        }


        Eigen::Array<Scalar, hmax, 1> vDiff;
        auto quadIt = quadData.begin();

        for (auto sector1d : matslise->sectors) {

            if (sector1d->backward) {
                vDiff = Eigen::Array<Scalar, hmax, 1>::Map(sector1d->vs.data());
                for (int i = 0; i < hmax; i += 2)
                    vDiff(i) *= -1;
            } else {
                vDiff = -Eigen::Array<Scalar, hmax, 1>::Map(sector1d->vs.data());
            }

            vDiff += Eigen::Array<Scalar, hmax, 1>::Map(
                    legendre::getCoefficients<hmax, Scalar, Scalar>(
                            [this, &y](const Scalar &x) -> Scalar { return this->sector2d->se2d->potential(x, y); },
                            sector1d->min, sector1d->max
                    ).data());

            Scalar h = 1;
            for (int i = 0; i < hmax; ++i, h *= sector1d->h)
                vDiff(i) *= h;

            for (int i = 0; i < N; ++i)
                for (int j = 0; j <= i; ++j) {
                    if (skipOdd && (i & 1) != (j & 1))
                        continue;
                    Scalar v = ((*quadIt) * vDiff).sum();
                    dV(i, j) += v;
                    if (i != j)
                        dV(j, i) += v;
                    ++quadIt;
                }
        }
        return dV;
    }

private:
    void calculateQuadData() {
#define hmax2 (2*MATSLISE_HMAX_delta-1)
        for (auto sector1d : matslise->sectors) {
            Eigen::Index N = sector2d->se2d->N;
            Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> u = sector1d->t_coeff.unaryExpr(
                    [](const Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign> &T) {
                        return T(0, 0);
                    });
            Eigen::Array<Scalar, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> v = sector1d->t_coeff.unaryExpr(
                    [](const Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign> &T) {
                        return T(0, 1);
                    });

            Eigen::Array<Eigen::Array<Scalar, hmax2, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> uu
                    = etaProduct<Scalar, hmax2>(u, u);
            Eigen::Array<Eigen::Array<Scalar, hmax2, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> uv
                    = etaProduct<Scalar, hmax2>(u, v);
            Eigen::Array<Eigen::Array<Scalar, hmax2, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta> vv
                    = etaProduct<Scalar, hmax2>(v, v);

            std::vector<matslise::Y<Scalar>> y0;
            for (auto &f : sector2d->eigenfunctions) {
                if (sector1d->backward) {
                    y0.emplace_back(f(sector1d->max));
                    y0.back().reverse();
                } else {
                    y0.emplace_back(f(sector1d->min));
                }
            }

            typedef Array<Scalar, Dynamic, 1> ArrayXs;
            ArrayXs grid = lobatto::grid<Scalar>(ArrayXs::LinSpaced(20, sector1d->min, sector1d->max));
            for (int i = 0; i < N; ++i)
                for (int j = 0; j <= i; ++j) {
                    if (skipOdd && (i & 1) != (j & 1))
                        continue;
                    Eigen::Array<Eigen::Array<Scalar, MATSLISE_INTEGRATE_delta, 1>, MATSLISE_ETA_delta, MATSLISE_ETA_delta>
                            eta = (i == j
                                   ? eta_integrals<Scalar, true>(sector1d->h,
                                                                 sector1d->vs[0] - sector2d->eigenvalues[i],
                                                                 sector1d->vs[0] - sector2d->eigenvalues[j])
                                   : eta_integrals<Scalar, false>(sector1d->h,
                                                                  sector1d->vs[0] - sector2d->eigenvalues[i],
                                                                  sector1d->vs[0] - sector2d->eigenvalues[j]));

                    /*
                    if (i == 9 && j == 9) {
                        std::cout << "\n****\n" << sector1d->h << ", "
                                  << sector1d->vs[0] - sector2d->eigenvalues[i] << ", "
                                  << sector1d->vs[0] - sector2d->eigenvalues[j] << std::endl;

                        Eigen::Array<Scalar, MATSLISE_ETA_delta * MATSLISE_ETA_delta, MATSLISE_INTEGRATE_delta> coutI;
                        for (int i = 0; i < MATSLISE_ETA_delta; ++i)
                            for (int j = 0; j < MATSLISE_ETA_delta; ++j)
                                coutI.row(i * MATSLISE_ETA_delta + j) = eta(i, j);
                        std::cout << coutI << std::endl;
                        std::cout << std::endl;
                    }
                    */

                    Eigen::Array<Scalar, hmax, 1> quadratures = Eigen::Array<Scalar, hmax, 1>::Zero();
                    {
                        Scalar suu = y0[i].y(0) * y0[j].y(0);
                        Scalar suv = y0[i].y(0) * y0[j].y(1);
                        Scalar svu = y0[i].y(1) * y0[j].y(0);
                        Scalar svv = y0[i].y(1) * y0[j].y(1);
                        for (int l = 0; l < MATSLISE_ETA_delta; ++l)
                            for (int m = 0; m < MATSLISE_ETA_delta; ++m) {
                                Eigen::Array<Scalar, hmax2, 1> s =
                                        (suu * uu(l, m) + suv * uv(l, m) + svu * uv(m, l) +
                                         svv * vv(l, m)).template head<hmax2>();
                                for (int n = std::max(0, 2 * l - 1 + 2 * m - 1);
                                     n < hmax2 && n < MATSLISE_INTEGRATE_delta; ++n) {
                                    int count = std::min(hmax, MATSLISE_INTEGRATE_delta - n);
                                    quadratures.segment(0, count) += s(n) * eta(l, m).segment(n, count);
                                }
                            }
                    }
                    legendreTransform<Scalar>(sector1d->h, quadratures);
                    quadData.push_back(quadratures);

                    if (sector1d->backward) {
                        auto &quad = quadData.back();
                        for (Eigen::Index k = 1; k < hmax; k += 2)
                            quad(k) *= -1;
                    }

                    /*
                    auto &quad = quadData.back();
                    Array<Scalar, Dynamic, 1> fi = sector2d->func_eigenfunctions[i](grid).unaryExpr(
                            [](const Y<Scalar> &y) { return y.y[0]; });
                    Array<Scalar, Dynamic, 1> fj = sector2d->func_eigenfunctions[j](grid).unaryExpr(
                            [](const Y<Scalar> &y) { return y.y[0]; });

                    if (abs(lobatto::quadrature<Scalar>(grid, fi * fj) - quad(0)) > 1e-3) {
                        cout << "u:\n" << u << endl;
                        cout << "\n *** quad <-> lobatto" << endl;
                        cout << "\n" << quad.transpose() << endl;
                        cout << "\n" << lobatto::quadrature<Scalar>(grid, fi * fj) << endl;
                        cout << "\n"
                             << lobatto::quadrature<Scalar>(grid,
                                                            fi * fj * (2 * (grid - sector1d->min) / sector1d->h - 1))
                             << endl;
                        cout << "\n"
                             << lobatto::quadrature<Scalar>(grid, fi * fj * (6 * (grid - sector1d->min) / sector1d->h *
                                                                             (grid - sector1d->min) / sector1d->h -
                                                                             6 * (grid - sector1d->min) / sector1d->h +
                                                                             1))
                             << endl;
                        cout << "wrong" << endl;
                    }
                    */

                }
        }
    }
};

template<typename Scalar, int hmax>
class BasisQuadratureHalf : public AbstractBasisQuadrature<Scalar> {
public:
    BasisQuadrature<Scalar, hmax, true> half;

    BasisQuadratureHalf(
            const typename matslise::Matslise2D<Scalar>::Sector *sector2d,
            const matslise::MatsliseHalf<Scalar> *matslise) {
        half = BasisQuadrature<Scalar, hmax, true>(sector2d, matslise->ms);
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dV(const Scalar &y) override {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dV = half.dV(y);
        dV *= 2;

        Eigen::Index N = half.sector2d->se2d->N;
        for (int i = 0; i < N; ++i)
            dV(i, i) -= half.sector2d->eigenvalues[i];

        return dV;
    }
};

#endif //MATSLISE_BASISQUADRATURE_H
