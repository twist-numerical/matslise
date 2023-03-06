#include <cmath>
#include "../test.h"
#include "../../matslise/liouville.h"
#include "../../matslise/util/constants.h"


using namespace matslise;
using namespace Eigen;

inline double square(double x) {
    return x * x;
}

TEST_CASE("Symbolic test liouville transformation", "[liouville]") {
    double rMin = 0.1;
    double rMax = 1.3;
    LiouvilleTransformation<double> transformation{
            {rMin, rMax},
            [](double r) { return std::cos(r); },
            [](double) { return 0; },
            [](double r) { return std::tan(r) * std::sin(r); }
    };

    auto dom = transformation.xDomain();
    REQUIRE_THAT(dom.min(), WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(dom.max(), WithinAbs(1.3136317360059335, 1e-6));

    for (double r = rMin; r < rMax; r += 0.01) {
        double exact = log(1 + square(std::tan(r))) / 2 - 0.005008355623235258;

        REQUIRE_THAT(transformation.r2x(r), WithinAbs(exact, 1e-5));
        REQUIRE_THAT(transformation.x2r(exact), WithinAbs(r, 1e-5));
    }

    for (double x = 0.1; x < dom.max(); x += 0.01) {
        double e2x = std::exp(2 * x);
        double e4x = e2x * e2x;
        double exactV = (-1.010067046422495 * e2x + 0.25) / (1.020235438268662 * e4x - 2.02013409284499 * e2x + 1);

        REQUIRE_THAT(transformation.V(x), WithinAbs(exactV, 1e-3));
    }
}
