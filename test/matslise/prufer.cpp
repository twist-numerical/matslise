#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

double zero(double x) {
    (void) x;
    return 0;
}

TEST_CASE("Prufer", "[matslise][prufer]") {
    Matslise<double> ms(&zero, 0, M_PI, 60);
    Y<double, 1> y0({1, -1}, {0, 0});
    REQUIRE(Approx(ms.propagate(16, y0, 0, M_PI / 2).second / M_PI -
                   ms.propagate(16, y0, M_PI, M_PI / 2).second / M_PI).margin(1e-10) == 4);
}
