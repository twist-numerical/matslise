#include "../catch.hpp"
#include "../../matslise/util/legendre.h"
#include "../../matslise/util/horner.h"

using namespace matslise::legendre;


TEST_CASE("As polynomial", "[util][legendre]") {
    double xmin = -0.5;
    double xmax = 2;

    Legendre<10, double> l{[](double x) { return std::sin(x); }, xmin, xmax};
    auto poly = l.asPolynomial();

    for (double x = xmin; x < xmax; x += 0.01) {
        REQUIRE(Approx(horner<double>(poly, (x - xmin) / (xmax - xmin), poly.size())).margin(1e-8) == std::sin(x));
    }
}
