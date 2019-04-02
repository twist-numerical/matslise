#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include <matslise/matslise.h>


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("V(x) = 0", "[matslise][simple]") {
    Matslise ms([](double x) -> double { return 0; }, -M_PI / 2, M_PI / 2, 11);


    Y<double> y0({0, 1}, {0, 0});

    for (int a = 1; a < 10; ++a) {
        double E = a * a;
        auto f = ms.eigenfunctionCalculator(E, y0, y0);
        double scale = (a % 2 == 1 ? f->eval(0) : f->eval(M_PI_2 / a)).y[0];
        for (double x = -M_PI_2 + 0.001; x < M_PI_2; x += 0.01) {
            REQUIRE(Approx(f->eval(x).y[0] / scale).margin(1e-7) == (a % 2 == 1 ? cos(a * x) : sin(a * x)));
        }
        delete f;
    }
}