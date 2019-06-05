#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "checkOrthonormality.h"
#include "../../src/se2d.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

void testQuartic(double a, double c, double alpha, double tolerance=0.0001) {
    SEnD<2> p(
            [a, c](double x, double y) -> double {
                return x * x + y * y + c * (x * x * x * x + 2 * a * x * x * y * y + y * y * y * y);
            },
            {{-alpha, alpha}, -alpha, alpha},
            Options<2>().tolerance(tolerance).stepsPerSector(2).N(10).nested(Options<1>().tolerance(tolerance)));
    p.findEigenvalue(6.5);
}

TEST_CASE("Eigenfunctions quartic: c=-3", "[se2d][eigenvalues][quartic]") {
    testQuartic(-1, 1e-3, 8);
    testQuartic(0, 1e-3, 8);
    testQuartic(1, 1e-3, 8);
}

TEST_CASE("Eigenfunctions quartic: c=3", "[se2d][eigenvalues][quartic]") {
    //testQuartic(-1, 1e3, 8);
    testQuartic(0, 1, 5, 0.0001);
    testQuartic(1, 1, 5, 0.0001);
}