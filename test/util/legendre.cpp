#include "../catch.hpp"
#include "../../matslise/util/legendre.h"

using namespace matslise::legendre;

TEST_CASE("As polynomial", "[util][legendre]") {
    double xmin = -0.5;
    double xmax = 2;

    Legendre<double, 10> l{[](double x) { return std::sin(x); }, xmin, xmax};
    auto poly = l.asPolynomial();

    for (double x = xmin; x < xmax; x += 0.01) {
        REQUIRE(Approx(poly((x - xmin) / (xmax - xmin))).margin(1e-8) == std::sin(x));
    }
}
