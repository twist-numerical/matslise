#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../matslise/matslise2d.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

const vector<tuple<Index, double, Index>> TWO_SEXTIC_EIGENVALUES{
        {0,  2 * 1.9922357643,  1},
        {1,  2 * 4.3051384550,  1},
        {2,  2 * 4.6993231357,  1},
        {3,  2 * 6.8954263765,  1},
        {4,  2 * 7.8378702941,  1},
        {5,  2 * 7.9593012390,  1},
        {6,  2 * 10.0165291976, 1},
        {7,  2 * 10.5861882834, 1},
        {8,  2 * 11.7788803250, 1},
        {9,  2 * 11.8005553313, 1},
        {10, 2 * 13.4155400229, 1},
        {11, 2 * 14.2097757808, 1},
        {12, 2 * 14.4819638906, 1}
};

double V6(double x) {
    return x * x * (.5 + x * x * (2 + .5 * x * x));
}

/* A hard problem
// http://www.sciencedirect.com/science/article/pii/S0021999196901400
TEST_CASE("Coupled anharmonic sextic oscillators", "[matslise2d][eigenfunctions][ixaru][auto]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-7;
    config.stepsPerSector = 3;
    config.ySectorBuilder = sector_builder::uniform<Matslise2D<>>(31);

    Matslise2D<> problem(
            [](double x, double y) -> double {
                return 2 * (V6(x) + V6(y) + x * y);
            }, {-4., 4., -4., 4.}, config);

    checkProblem(problem, TWO_SEXTIC_EIGENVALUES);
}
 */