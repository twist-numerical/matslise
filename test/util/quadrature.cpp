#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../matslise/matslise.h"
#include "../../matslise/util/quadrature.h"


using namespace std;
using namespace Eigen;
using namespace quadrature;


int max(int a, int b) {
    return a > b ? a : b;
}

void testFunction(const function<double(double)> &f, const function<double(double)> &integral, double a, double b, int Nmin = 10) {
    for (int N = Nmin; N < 3 * Nmin || N < 100; N += 10) {
        ArrayXd grid(N);
        for (int i = 0; i < N; ++i)
            grid[i] = a + (b - a) * i / (N - 1);
        ArrayXd lg = lobatto::grid(grid);
        ArrayXd eval = lg.unaryExpr(f);
        CHECK(Approx(lobatto::quadrature<double>(lg, eval)).margin(1e-10) == (integral(b) - integral(a)));
    }
}

void testGaussKronrod(const function<double(double)> &f, const function<double(double)> &integral, double a, double b) {
    CHECK(Approx(gauss_kronrod::adaptive(f, a, b, 1e-10)).margin(1e-10) == (integral(b) - integral(a)));
}

TEST_CASE("adaptive sin(x) (gaussKronrod)", "[util][gaussKronrod]") {
    for (double a = -20; a < 10; a += 3.94) {
        for (double b = a + .1; b < 15; b += 3.941) {
            testGaussKronrod(
                    [](double x) { return sin(x); },
                    [](double x) { return -cos(x); },
                    a, b);
        }
    }
}

TEST_CASE("sin(x)", "[util][lobatto]") {
    for (double a = -20; a < 10; a += 3.94) {
        for (double b = a + .1; b < 15; b += 3.941) {
            testFunction(
                    [](double x) { return sin(x); },
                    [](double x) { return -cos(x); },
                    a, b, max(static_cast<int>(b - a), 10));
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
