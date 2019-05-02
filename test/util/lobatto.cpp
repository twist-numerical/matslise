#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matscs.h"
#include "../../src/util/lobatto.h"


using namespace std;
using namespace Eigen;


int max(int a, int b) {
    return a > b ? a : b;
}

void testFunction(function<double(double)> f, function<double(double)> integral, double a, double b, int Nmin = 10) {
    for (int N = Nmin; N < 3 * Nmin || N < 100; N += 10) {
        ArrayXd grid(N);
        for (int i = 0; i < N; ++i)
            grid[i] = a + (b - a) * i / (N - 1);
        ArrayXd lg = lobatto::grid(grid);
        ArrayXd eval = lobatto::apply<1>(&lg, f);
        CHECK(Approx(lobatto::multi_quadrature<1>(&lg, eval)).margin(1e-10) == (integral(b) - integral(a)));
    }
}

TEST_CASE("sin(x)", "[util][lobatto]") {
    for (double a = -20; a < 10; a += 3.94) {
        for (double b = a + .1; b < 15; b += 3.941) {
            testFunction(
                    [](double x) { return sin(x); },
                    [](double x) { return -cos(x); },
                    a, b, max((int) b - a, 10));
        }
    }
}

TEST_CASE("sin(c*x)", "[util][lobatto]") {
    for (double a = -10; a < 10; a += 3.94) {
        for (double b = a + .1; b < 15; b += 3.941) {
            for (double c = .002; c < 6; c *= 1.12)
                testFunction(
                        [c](double x) { return sin(c * x); },
                        [c](double x) { return -cos(c * x) / c; },
                        a, b, max(20, (int) (c * 30)));
        }
    }
}


TEST_CASE("cos(x)sin(c*x)", "[util][lobatto]") {
    for (double a = -6; a < 10; a += 3.94) {
        for (double b = a + .1; b < 11; b += 3.941) {
            for (double c = .02; c < 5; c *= 1.12)
                testFunction(
                        [c](double x) { return cos(x) * sin(c * x); },
                        [c](double x) { return -0.5 * cos((c + 1) * x) / (c + 1) - 0.5 * cos((c - 1) * x) / (c - 1); },
                        a, b, max(30, (int) (c * 30)));
        }
    }
}

