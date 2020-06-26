#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Eigenfunctions henon", "[matslise2d][eigenfunctions][henon]") {
    Matslise2D<> p2(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * x * (y * y - x * x / 3);
            },
            {{-6, 6}, -6, 6},
            Options2<>().sectorCount(21).nested(Options1<>().sectorCount(16)));
    pair<double, int> eigenvalues[] = {
            {0.998594690530479, 1},
            {1.99007660445524,  2},
            {2.95624333869018,  1},
            {2.98532593386986,  2},
            {3.92596412795287,  2},
            {3.98241882458866,  1},
            {3.98575763690663,  1},
            {4.87014557482289,  1},
            {4.89864497284387,  2}};

    CHECK(Approx(p2.firstEigenvalue()).margin(1e-7) == 2 * eigenvalues[0].first);

    int n = 2;
    ArrayXd x = ArrayXd::LinSpaced(n, -5, 5);
    double E;
    int multiplicity;
    for (auto &Em  : eigenvalues) {
        tie(E, multiplicity) = Em;
        E *= 2;

        const pair<double, double> &error = p2.matchingError(E);
        REQUIRE(abs(error.first) < 1e-3);
        REQUIRE(Approx(p2.eigenvalue(E)).margin(1e-7) == E);

        const vector<Eigenfunction2D<>> f = p2.eigenfunction(E);
        REQUIRE(f.size() == multiplicity);

        for (int k = 0; k < multiplicity; ++k) {
            auto fkxx = f[k](x, x);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    REQUIRE(Approx(f[k](x[i], x[j])).margin(1e-7) == fkxx(i, j));
        }
    }
}

TEST_CASE("Eigenfunctions henon (half)", "[matslise2d][eigenfunctions][henon][half]") {
    Matslise2D<> p2(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * y * (x * x - y * y / 3);
            },
            {{-6, 6}, -6, 6},
            Options2<>().sectorCount(21).nested(Options1<>().symmetric(true).sectorCount(8)));
    pair<double, int> eigenvalues[] = {
            {0.998594690530479, 1},
            {1.99007660445524,  2},
            {2.95624333869018,  1},
            {2.98532593386986,  2},
            {3.92596412795287,  2},
            {3.98241882458866,  1},
            {3.98575763690663,  1},
            {4.87014557482289,  1},
            {4.89864497284387,  2}};

    CHECK(Approx(p2.firstEigenvalue()).margin(1e-7) == 2 * eigenvalues[0].first);

    int n = 2;
    ArrayXd x = ArrayXd::LinSpaced(n, -5, 5);
    double E;
    int multiplicity;
    for (auto &Em  : eigenvalues) {
        tie(E, multiplicity) = Em;
        E *= 2;

        const pair<double, double> &error = p2.matchingError(E);
        REQUIRE(abs(error.first) < 1e-3);
        REQUIRE(Approx(p2.eigenvalue(E)).margin(1e-7) == E);

        const vector<Eigenfunction2D<>> f = p2.eigenfunction(E);
        REQUIRE(f.size() == multiplicity);

        for (int k = 0; k < multiplicity; ++k) {
            auto fkxx = f[k](x, x);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    REQUIRE(Approx(f[k](x[i], x[j])).margin(1e-7) == fkxx(i, j));
        }
    }
}
