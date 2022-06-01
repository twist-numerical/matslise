#include "../catch.hpp"
#include "../../matslise/util/legendre.h"
#include "fEquals.h"

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

TEST_CASE("As polynomial with derivatives", "[util][legendre]") {
    double xmin = -0.1;
    double xmax = 1.5;
    double h = xmax - xmin;

    Legendre<double, 14> l{[](double x) { return std::cos(x); }, xmin, xmax};

    auto scaled = [&](const auto &f) {
        return [&, f](double x) { return f(xmin + h * x); };
    };

    fEquals<double>(scaled([](double x) { return std::cos(x); }),
                    l.asPolynomial(), 0, 1);
    fEquals<double>(scaled([](double x) { return -std::sin(x); }),
                    l.asPolynomial().derivative() / h, 0, 1);
    fEquals<double>(scaled([](double x) { return -std::cos(x); }),
                    l.asPolynomial().template derivative<2>() / (h * h), 0, 1);
    fEquals<double>(scaled([](double x) { return std::sin(x); }),
                    l.asPolynomial().template derivative<3>() / (h * h * h), 0, 1);
    fEquals<double>(scaled([](double x) { return std::cos(x); }),
                    l.asPolynomial().template derivative<4>() / (h * h * h * h), 0, 1);
}
