#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("V(x) = 0", "[matslise][simple]") {
    Matslise<double> ms([](double x) -> double {
        (void) x;
        return 0;
    }, -constants<double>::PI / 2, constants<double>::PI / 2, 11);


    Y<double> y0({0, 1}, {0, 0});

    for (int a = 1; a < 10; ++a) {
        double E = a * a;
        std::function<Y<double>(double)> f = ms.eigenfunctionCalculator(E, y0, y0);
        double scale = (a % 2 == 1 ? f(0) : f(constants<double>::PI / 2 / a)).y[0];
        for (double x = -constants<double>::PI / 2 + 0.001; x < constants<double>::PI / 2; x += 0.01) {
            REQUIRE(Approx(f(x).y[0] / scale).margin(1e-7) == (a % 2 == 1 ? cos(a * x) : sin(a * x)));
        }
    }
}

TEST_CASE("V(x) = 0 (auto)", "[matslise][simple][auto)") {
    Matslise<double> ms([](double x) -> double {
        (void) x;
        return 0;
    }, -constants<double>::PI / 2, constants<double>::PI / 2, Matslise<double>::AUTO(1e-8));


    Y<double> y0({0, 1}, {0, 0});

    for (int a = 1; a < 10; ++a) {
        double E = a * a;
        std::function<Y<double>(double)> f = ms.eigenfunctionCalculator(E, y0, y0);
        double scale = (a % 2 == 1 ? f(0) : f(constants<double>::PI / 2 / a)).y[0];
        for (double x = -constants<double>::PI / 2 + 0.001; x < constants<double>::PI / 2; x += 0.01) {
            REQUIRE(Approx(f(x).y[0] / scale).margin(1e-7) == (a % 2 == 1 ? cos(a * x) : sin(a * x)));
        }
    }
}