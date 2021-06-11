#include <cmath>
#include <set>
#include <tuple>
#include "../catch.hpp"
#include "../../matslise/matslise2d.h"
#include "../../matslise/util/quadrature.h"
#include "checkOrthonormality.h"
#include "../../matslise/util/constants.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Eigenvalues V=0", "[matslise2d][eigenfunctions][zero]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-7;

    Matslise2D<> p(
            [](double, double) -> double {
                return 0;
            },
            {0., constants<double>::PI, 0., constants<double>::PI}, config);

    checkProblem(p, vector<tuple<Index, double, Index>>{
            {0,  2,  1},
            {1,  5,  2},
            {3,  8,  1},
            {4,  10, 2},
            {6,  13, 2},
            {8,  17, 2},
            {10, 18, 1},
            {11, 20, 2},
            {13, 25, 2},
            {15, 26, 2},
            {17, 29, 2},
            {19, 32, 1},
    });
}
