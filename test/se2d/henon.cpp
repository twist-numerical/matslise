#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "./checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

const vector<tuple<Index, double, Index>> HENON_EIGENVALUES{
        {0,  2 * 0.998594690530479, 1},
        {1,  2 * 1.99007660445524,  2},
        {3,  2 * 2.95624333869018,  1},
        {4,  2 * 2.98532593386986,  2},
        {6,  2 * 3.92596412795287,  2},
        {8,  2 * 3.98241882458866,  1},
        {9,  2 * 3.98575763690663,  1},
        {10, 2 * 4.87014557482289,  1},
        {11, 2 * 4.89864497284387,  2}};

TEST_CASE("Eigenfunctions henon", "[matslise2d][eigenfunctions][henon][auto]") {
    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * x * (y * y - x * x / 3);
            },
            {{-6, 6}, -6, 6},
            Options2<>().tolerance(1e-6).nested(Options1<>().tolerance(1e-7)));

    checkProblem(problem, HENON_EIGENVALUES);
}

TEST_CASE("Eigenfunctions henon (symmetric)", "[matslise2d][eigenfunctions][henon][symmetric][auto]") {
    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * y * (x * x - y * y / 3);
            },
            {{-6, 6}, -6, 6},
            Options2<>().tolerance(1e-6).nested(Options1<>().tolerance(1e-7)));

    checkProblem(problem, HENON_EIGENVALUES);
}

TEST_CASE("Eigenfunctions henon (half)", "[matslise2d][eigenfunctions][henon][half][auto]") {
    Matslise2DHalf<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * x * (y * y - x * x / 3);
            },
            {{-6, 6}, -6, 6},
            Options2<>().tolerance(1e-6).nested(Options1<>().tolerance(1e-7)));

    checkProblem(problem, HENON_EIGENVALUES);
}
