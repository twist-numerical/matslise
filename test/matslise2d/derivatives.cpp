#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../matslise/matslise.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Value>
Value approxDiff(const std::function<Value(double)> &f, double x, double h = 1e-5) {
    return (f(x + h) - f(x - h)) / (2 * h);
}

TEST_CASE("E derivative of matchingErrorMatrix", "[matslise2d][derivatives][e_derivative]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-6;

    Matslise2D<> p(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {-5., 5., -5., 5.}, config);

    MatrixXd mat, diff;
    for (double E = 0.1; E < 20; E += 0.2) {
        INFO("E: " << E);
        tie(mat, diff) = p.matchingErrorMatrix(E);
        MatrixXd adiff = approxDiff<MatrixXd>([&](double E) -> MatrixXd { return p.matchingErrorMatrix(E).first; }, E,
                                              1e-6);
        CHECK(((diff - adiff).array().abs() / (diff.array().abs().max(1.))).maxCoeff() < 1e-4);
    }
}

TEST_CASE("Derivatives of eigenfunctions", "[matslise2d][y_derivative][derivatives][x_derivative]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-6;

    Matslise2D<> p(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {-5., 5., -5., 5.}, config);

    double step = .4;
    double phi, phi_x, phi_y, a_phi_x, a_phi_y;
    for (double E : {3.19, 5.5, 7.55, 8., 8.45, 9.9, 11.3, 12.1, 12.2, 13.9}) {
        E = p.eigenvalue(E).first;
        INFO("E: " << E);
        for (const auto &f : p.eigenfunctionWithDerivatives(E))
            for (double x = -5 + step / 2; x < 5; x += step)
                for (double y = -5 + step / 2; y < 5; y += step) {
                    INFO("(x, y): (" << x << ", " << y << ")");
                    tie(phi, phi_x, phi_y) = f(x, y);
                    a_phi_x = approxDiff<double>([&](double x) -> double {
                        return get<0>(f(x, y));
                    }, x, 1e-6);
                    a_phi_y = approxDiff<double>([&](double y) -> double {
                        return get<0>(f(x, y));
                    }, y, 1e-7);

                    CHECK(Approx(phi_x).margin(1e-4) == a_phi_x);
                    CHECK(Approx(phi_y).margin(1e-4) == a_phi_y);
                }
    }
}
