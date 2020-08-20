#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/util/constants.h"
#include "../../src/matslise.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

// Numerical approximation (estimated error < 2e-7)
const vector<tuple<Index, double, Index>> PENTAGON_EIGENVALUES{
        {0,  2.023972511137453,  1},
        // Finding orthogonal eigenfunctions is hard
        {1,  4.423508358409953,  2},
        {3,  7.01457283417795,   1},
        {4,  7.019042961355106,  1},
        {5,  7.305044641507928,  1},
        {6,  9.89999119875196,   1},
        {7,  9.90362055869358,   1},
        {8,  10.149917029739067, 1},
        {9,  10.151510900699652, 1},
        {10, 12.937852419345413, 1}
};

TEST_CASE("2D: Pentagon", "[matslise2d][eigenfunctions][slow]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-7;
    config.basisSize = 14;
    config.ySectorBuilder = sector_builder::automatic<Matslise2D<>, true>(1e-6);

    Matslise2D<> problem(
            [](double x, double y) -> double {
                const double pi = constants<double>::PI;
                double arg = atan2(y, x);
                if (arg < 0)
                    arg += 2. * pi;
                double r = 1. / cos((arg - floor(arg / (pi * 2. / 5.)) * (pi * 2. / 5.) - pi / 5.));
                return pow(sqrt(x * x + y * y) / r, 3.);
            },
            {-3., 3., -3., 3.}, config);

    checkProblem(problem, PENTAGON_EIGENVALUES, 1e-3);
}
