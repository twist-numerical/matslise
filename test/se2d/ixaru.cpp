#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

const vector<tuple<Index, double, Index>> IXARU_EIGENVALUES{
        {0,  3.1959181,  1},
        {1,  5.5267439,  2},
        {3,  7.5578033,  1},
        {4,  8.0312723,  1},
        {5,  8.4445814,  1},
        {6,  9.9280611,  2},
        {8,  11.3118171, 2},
        {10, 12.1032536, 1},
        {11, 12.2011790, 1},
        {12, 13.3323313, 1}
};

TEST_CASE("Eigenfunctions ixaru auto", "[matslise2d][eigenfunctions][ixaru][auto]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-6;

    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {-5.5, 5.5, -5.5, 5.5}, config);

    checkProblem(problem, IXARU_EIGENVALUES);
}

TEST_CASE("Eigenfunctions ixaru halfrange", "[matslise2d][eigenfunctions][ixaru][halfrange]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-6;

    Matslise2DHalf<> problem(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {-5.5, 5.5, -5.5, 5.5}, config);

    checkProblem(problem, IXARU_EIGENVALUES);
}

TEST_CASE("Eigenfunctions ixaru auto high n", "[matslise2d][eigenfunctions][ixaru][auto][slow]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-6;
    config.basisSize = 20;

    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {-5.5, 5.5, -5.5, 5.5}, config);

    checkProblem(problem, IXARU_EIGENVALUES);
}
