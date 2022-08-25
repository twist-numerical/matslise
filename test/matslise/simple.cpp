#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../matslise/matslise.h"
#include "../../matslise/util/constants.h"


using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace matslise::sector_builder;

TEST_CASE("V(x) = 0", "[matslise][simple]") {
    Matslise<double> ms([](double x) -> double {
        (void) x;
        return 0;
    }, -constants<double>::PI / 2, constants<double>::PI / 2, 1e-8, UniformSectorBuilder<Matslise<>>(11));


    Y<double> y0({0, 1}, {0, 0});

    for (int a = 1; a < 50; ++a) {
        double E = a * a;
        unique_ptr<Matslise<>::Eigenfunction> f = ms.eigenfunction(E, y0, y0);
        double scale = (a % 2 == 1 ? (*f)(0) : (*f)(constants<double>::PI / 2 / a))[0];
        for (double x = -constants<double>::PI / 2 + 0.001; x < constants<double>::PI / 2; x += 0.01) {
            REQUIRE(Approx(abs((*f)(x)[0])).margin(1e-7) == abs(scale * (a % 2 == 1 ? cos(a * x) : sin(a * x))));
        }
    }
}

TEST_CASE("V(x) = 0 (auto)", "[matslise][simple][auto]") {
    Matslise<double> ms([](double x) -> double {
        (void) x;
        return 0;
    }, 0, constants<double>::PI, 1e-6);


    Y<double> y0({0, 1}, {0, 0});

    for (int a = 1; a < 50; ++a) {
        double E = a * a;
        REQUIRE(ms.eigenvalueError(E, y0, y0) < 1e-6);
        unique_ptr<Matslise<>::Eigenfunction> f = ms.eigenfunction(E, y0, y0);
        double scale = sqrt(2 / constants<double>::PI);
        for (double x = 0.005; x < constants<double>::PI; x += 0.01) {
            REQUIRE(Approx(abs((*f)(x)[0])).margin(1e-7) == abs(scale * sin(a * x)));
        }
    }
}