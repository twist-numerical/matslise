#include <cmath>
#include <vector>
#include <tuple>
#include "../test.h"
#include "../../matslise/matslise.h"
#include "../../matslise/util/constants.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

double zero(double x) {
    (void) x;
    return 0;
}

TEST_CASE("Prufer", "[matslise][prufer]") {
    Matslise<double> ms(&zero, 0, constants<double>::PI, 1e-8);
    Y<double, 1> y0({1, -1}, {0, 0});
    REQUIRE_THAT(ms.propagate(16, y0, 0, constants<double>::PI / 2).second / constants<double>::PI -
                 ms.propagate(16, y0, constants<double>::PI, constants<double>::PI / 2).second /
                 constants<double>::PI, WithinAbs(4.0, 1e-10));
}
