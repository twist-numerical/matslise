#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../matslise/util/legendre.h"
#include "../../matslise/matslise2d.h"
#include "../../matslise/matslise2d/etaIntegrals.h"
#include "../../matslise/matslise2d/basisQuadrature.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
typename Matslise2D<Scalar>::MatrixXs
calculateDeltaV(const typename Matslise2D<Scalar>::Sector *sector, const Scalar &y) {
    const Matslise2D<Scalar> *se2d = sector->se2d;
    Index N = se2d->config.basisSize;
    typename Matslise2D<Scalar>::MatrixXs dV(N, N);

    const Scalar ybar = sector->ybar;
    function<Scalar(Scalar)> vbar_fun = [se2d, ybar](Scalar x) -> Scalar { return se2d->potential(x, ybar); };

    typename Matslise2D<Scalar>::ArrayXs grid = lobatto::grid<Scalar>(
            Matslise2D<Scalar>::ArrayXs::LinSpaced(511, se2d->domain.min(0), se2d->domain.max(0)));
    typename Matslise2D<Scalar>::ArrayXs vbar = grid.unaryExpr(vbar_fun);

    typename Matslise2D<Scalar>::ArrayXs vDiff =
            grid.unaryExpr([sector, y](const Scalar &x) -> Scalar { return sector->se2d->potential(x, y); }) -
            vbar;

    vector<typename Matslise2D<Scalar>::ArrayXs> eigenfunctions;
    eigenfunctions.reserve(N);

    for (int i = 0; i < N; ++i) {
        eigenfunctions.push_back(
                sector->eigenfunctions[i](grid).unaryExpr([](const Y<Scalar> &y) { return y.data(0); }));
        for (int j = 0; j <= i; ++j) {
            dV(i, j) = lobatto::quadrature<Scalar>(
                    grid,
                    eigenfunctions[i]
                    * vDiff * eigenfunctions[j]);
            if (i == j) {
                Scalar n = lobatto::quadrature<Scalar>(grid, eigenfunctions[i] * eigenfunctions[j]);
                if (abs(n - 1) > 1e-4)
                    cout << "Not normalized!! " << n << endl;
            }
            if (j < i) dV(j, i) = dV(i, j);
        }
        dV(i, i) += sector->eigenvalues[i];
    }

    return dV;
}

double Q_V6(double x) {
    return x * x * (.5 + x * x * (2 + .5 * x * x));
}

double Q_V(double x, double y) {
    return 2 * (Q_V6(x) + Q_V6(y) + x * y);
}

/* Some numerical nonsense when |v_1| >> 0
template<bool equal>
void checkEtaIntegral(double h, double dZ1, double dZ2) {
    Array<double, MATSLISE_INTEGRATE_eta_rows, Eigen::Dynamic> integrals
            = eta_integrals<double, equal>(h, dZ1, dZ2);

    ArrayXd xs = lobatto::grid<double>(ArrayXd::LinSpaced(1021, 0, h));

    Array<double, MATSLISE_ETA_delta, Dynamic> eta_i = calculateEta<double, MATSLISE_ETA_delta>(dZ1 * xs * xs);
    Array<double, MATSLISE_ETA_delta, Dynamic> eta_j = calculateEta<double, MATSLISE_ETA_delta>(dZ2 * xs * xs);

    Array<int, MATSLISE_ETA_delta, MATSLISE_ETA_delta> valid;
    valid.block(0, 0, 2, 2) << 0, 1, 1, 2;
    for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
        if (i > 1) {
            valid(0, i) = 2 + max(valid(0, i - 1), valid(0, i - 2));
            valid(1, i) = 2 + max(valid(1, i - 1), valid(1, i - 2));
        }
        for (int j = 2; j < MATSLISE_ETA_delta; ++j) {
            valid(j, i) = 2 + max(valid(j - 1, i), valid(j - 2, i));
        }
    }

    for (int i = 0; i < MATSLISE_ETA_delta; ++i)
        for (int j = 0; j < MATSLISE_ETA_delta; ++j) {
            Index ij = ETA_index(i, j);
            if (ij >= MATSLISE_INTEGRATE_eta_rows) break;
            ArrayXd xsn = xs.unaryExpr([](double) -> double { return 1; });
            for (int n = 0; n < 10; ++n) {
                if (n >= valid(i, j)) {
                    double quad = lobatto::quadrature<double>(xs, (eta_i.row(i) * eta_j.row(j)).transpose() * xsn);
                    double symb = integrals(ij, n);
                    INFO("h: " << h << ", dZ1: " << dZ1 << ", dZ2: " << dZ2 << ", i: " << i << ", j: " << j << ", n:"
                               << n << ", approx: " << quad << ", symbolic: " << symb);
                    CHECK((Approx(quad).epsilon(1e-2) == symb || Approx(quad).margin(1e-10) == symb));
                }
                xsn *= xs;
            }
        }
}


TEST_CASE("BasisQuadratures", "[matslise2d]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-7;
    config.stepsPerSector = 1;
    config.ySectorBuilder = sector_builder::uniform<Matslise2D<>>(31);

    Matslise2D<> p2d(
            [](double x, double y) -> double {
                return 2 * (Q_V6(x) + Q_V6(y) + x * y);
            }, {-4., 4., -4., 4.}, config);

    for (auto &sector : p2d.sectors) {
        const Index &N = config.basisSize;


        ArrayXd ys = ArrayXd::LinSpaced(5, sector->min + 1e-6, sector->max - 1e-6);
        for (Index yi = 0; yi < ys.rows(); ++yi) {
            double y = ys[yi];
            MatrixXd deltaV_approx = calculateDeltaV(sector, y);
            MatrixXd deltaV_symb = sector->quadratures->dV(*sector, y);

            if ((deltaV_approx - deltaV_symb).array().abs().maxCoeff() > 10) {
                cout << "y: " << y << " | yi: " << yi << endl;
                // cout << (deltaV_approx - deltaV_symb) << "\n" << endl;
                // cout << deltaV_approx << "\n" << endl;
                // cout << deltaV_symb << "\n\n" << endl;
            }
        }

        auto &quadData = dynamic_cast<BasisQuadrature<double, MATSLISE2D_DELTA_V_DEGREE> &>(*sector->quadratures).quadData;
        auto quadIterator = quadData.begin();
        for (auto &sector1d : dynamic_cast<Matslise<> &>(*sector->matslise).sectors) {
            ArrayXd xs = lobatto::grid<double>(ArrayXd::LinSpaced(101, sector1d->min, sector1d->max));
            vector<ArrayXd> eigenfunctions;
            for (int i = 0; i < N; ++i) {
                eigenfunctions.emplace_back(
                        sector->eigenfunctions[i](xs).unaryExpr([](const Y<> &y) { return y.y()[0]; }));
            }

            const int LEGENDRE_POLYNOMIALS_TO_CHECK = 1;
            ArrayXd xsz = (xs - sector1d->min) / (sector1d->max - sector1d->min);
            array<ArrayXd, 3> legendre_polynomials{
                    xs.unaryExpr([](double) -> double { return 1.; }),
                    2 * xsz - 1,
                    6 * xsz * xsz - 6 * xsz + 1
            };

            for (int i = 0; i < N; ++i)
                for (int j = 0; j <= i; ++j) {
                    Array<double, MATSLISE2D_DELTA_V_DEGREE, 1> &data = *quadIterator;

                    bool wrong = false;

                    for (int poly_i = 0; poly_i < LEGENDRE_POLYNOMIALS_TO_CHECK; ++poly_i) {
                        double quad = lobatto::quadrature<double>(xs, legendre_polynomials[poly_i] * eigenfunctions[i] *
                                                                      eigenfunctions[j]);
                        if (abs(quad - data[poly_i]) > 0.1) {
                            cout << "\n\n" << ArrayXd::Map(sector->eigenvalues.data(), N).transpose() << endl;
                            wrong = true;
                            cout << "poly: " << poly_i << " | i: " << i << " | j: " << j << " | xmin: " << sector1d->min
                                 << " | xmax: " << sector1d->max << endl;
                            cout << quad << " <-> " << data[poly_i] << "\n" << endl;
                            cout << "vs: " << ArrayXd::Map(sector1d->vs.data(), 10).transpose() << endl;

                            Array<double, MATSLISE_INTEGRATE_eta_rows, Eigen::Dynamic> integrals;

                            if (i == j) {
                                integrals = eta_integrals<double, true>(sector1d->h,
                                                                        sector1d->vs[0] - sector->eigenvalues[i],
                                                                        sector1d->vs[0] - sector->eigenvalues[j]);
                                checkEtaIntegral<true>(sector1d->h, sector1d->vs[0] - sector->eigenvalues[i],
                                                       sector1d->vs[0] - sector->eigenvalues[j]);
                            } else {
                                integrals = eta_integrals<double, false>(sector1d->h,
                                                                         sector1d->vs[0] - sector->eigenvalues[i],
                                                                         sector1d->vs[0] - sector->eigenvalues[j]);
                                checkEtaIntegral<false>(sector1d->h, sector1d->vs[0] - sector->eigenvalues[i],
                                                        sector1d->vs[0] - sector->eigenvalues[j]);
                            }


                            Array<double, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> u = sector1d->t_coeff.unaryExpr(
                                    [](const Matrix<double, 2, 2, DontAlign> &T) {
                                        return T(0, 0);
                                    });
                            Array<double, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> v = sector1d->t_coeff.unaryExpr(
                                    [](const Matrix<double, 2, 2, DontAlign> &T) {
                                        return T(0, 1);
                                    });

                            Array<double, MATSLISE_INTEGRATE_eta_rows, Dynamic>
                                    uu = etaProduct<double>(u, u);
                            Array<double, MATSLISE_INTEGRATE_eta_rows, Dynamic>
                                    uv = etaProduct<double>(u, v);
                            Array<double, MATSLISE_INTEGRATE_eta_rows, Dynamic>
                                    vv = etaProduct<double>(v, v);
                            Array<double, MATSLISE_INTEGRATE_eta_rows, Dynamic>
                                    vu(MATSLISE_INTEGRATE_eta_rows, MATSLISE_INTEGRATE_delta);

                            {
                                for (Index i = 0; i < MATSLISE_ETA_delta; ++i)
                                    for (Index j = 0; j < MATSLISE_ETA_delta; ++j) {
                                        Index ij = ETA_index(i, j);
                                        Index ji = ETA_index(j, i);
                                        if (ij >= MATSLISE_INTEGRATE_eta_rows ||
                                            ji >= MATSLISE_INTEGRATE_eta_rows)
                                            break;
                                        vu.row(ji) = uv.row(ij);
                                    }
                            }

                            cout << "\nvu:\n" << vu << endl;

                            for (Index i = 0; i < MATSLISE_ETA_delta; ++i)
                                for (Index j = 0; j < MATSLISE_ETA_delta; ++j) {
                                    Index ij = ETA_index(i, j);
                                    Index ji = ETA_index(j, i);
                                    if (ij < MATSLISE_INTEGRATE_eta_rows && ji < MATSLISE_INTEGRATE_eta_rows) break;
                                    vu.row(ji) = uv.row(ij);
                                }


                            Y<> yi = sector->eigenfunctions[i](
                                    sector1d->direction == Direction::forward ? sector1d->min : sector1d->max);
                            Y<> yj = sector->eigenfunctions[j](
                                    sector1d->direction == Direction::forward ? sector1d->min : sector1d->max);
                            if (sector1d->direction == Direction::backward) {
                                yi.reverse();
                                yj.reverse();
                            }

                            // cout << "Initial 0,0: " << sector->eigenfunctions[0](sector1d->direction == Direction::forward ? sector1d->min : sector1d->max).y().transpose() << endl;

                            double suu = yi.y()[0] * yj.y()[0];
                            double suv = yi.y()[0] * yj.y()[1];
                            double svu = yi.y()[1] * yj.y()[0];
                            double svv = yi.y()[1] * yj.y()[1];

                            cout << "direction: " << (sector->direction == Direction::forward ? "forward" : "backward")
                                 << endl;
                            cout << "initial: " << yi.y().transpose() << ", " << yj.y().transpose() << endl;
                            cout << "factors: " << suu << ", " << suv << ", " << svu << ", " << svv << endl;

                            Array<double, MATSLISE_INTEGRATE_eta_rows, Dynamic> fifj
                                    = (suu * uu + suv * uv + svu * vu + svv * vv);
                            Array<double, MATSLISE_INTEGRATE_eta_rows, Dynamic> int_fifj = integrals * fifj;

                            Array<int, MATSLISE_ETA_delta, MATSLISE_ETA_delta> valid;
                            valid.block(0, 0, 2, 2) << 0, 1, 1, 2;
                            for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
                                if (i > 1) {
                                    valid(0, i) = 2 + max(valid(0, i - 1), valid(0, i - 2));
                                    valid(1, i) = 2 + max(valid(1, i - 1), valid(1, i - 2));
                                }
                                for (int j = 2; j < MATSLISE_ETA_delta; ++j) {
                                    valid(j, i) = 2 + max(valid(j - 1, i), valid(j - 2, i));
                                }
                            }
                            for (int i = 0; i < MATSLISE_ETA_delta; ++i) {
                                for (int j = 0; j < MATSLISE_ETA_delta; ++j) {
                                    Index ij = ETA_index(i, j);
                                    if (ij >= MATSLISE_INTEGRATE_eta_rows) break;
                                    for (int n = 0; n < valid(i, j) && n < MATSLISE_INTEGRATE_delta; ++n)
                                        CHECK(int_fifj(ij, n) == 0);
                                }
                            }

                            //  cout << "integrals: (" << integrals.rows() << ", " << integrals.cols() << ")" << endl;
                            cout << "(integrals * fifj).sum(): " << int_fifj.sum() << endl;
                            cout << "(integrals * fifj).abs().maxCoeff(): " << int_fifj.abs().maxCoeff()
                                 << endl;
                            cout << "integrals * fifj:\n" << int_fifj << "\n" << endl;
                            cout << "fifj:\n" << fifj << "\n" << endl;
                            cout << "integrals:\n" << integrals << "\n" << endl;
                            cout << "u:\n" << u << "\n" << endl;
                            cout << "v:\n" << v << "\n" << endl;
                        }
                    }

                    if (wrong) {
                        cout << "h: " << sector1d->h << " | Zi: " << (sector1d->vs[0] - sector->eigenvalues[i])
                             << " | Zj: " << (sector1d->vs[0] - sector->eigenvalues[j]) << endl;
                        REQUIRE(false);
                    }

                    ++quadIterator;
                }

        }
    }
}


TEST_CASE("etaIntegrals", "[matslise2d]") {
    checkEtaIntegral<true>(0.821683, -0.791717, -0.791717);

    for (double h = 0.1; h < 1; h *= 1.3) {
        for (double Z1 = .01; Z1 < 50; Z1 *= 3) {
            checkEtaIntegral<true>(h, Z1, Z1);
            checkEtaIntegral<true>(h, -Z1, -Z1);

            for (double Z2 = Z1 * 1.4; Z2 < 50; Z2 *= 3) {
                checkEtaIntegral<false>(h, Z1, Z2);
                checkEtaIntegral<false>(h, Z2, Z1);
                checkEtaIntegral<false>(h, Z1, -Z2);
                checkEtaIntegral<false>(h, Z2, -Z1);
                checkEtaIntegral<false>(h, -Z1, Z2);
                checkEtaIntegral<false>(h, -Z2, Z1);
                checkEtaIntegral<false>(h, -Z1, -Z2);
                checkEtaIntegral<false>(h, -Z2, -Z1);
            }
        }
    }
}
 */