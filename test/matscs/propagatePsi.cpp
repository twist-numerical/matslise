#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

inline double square(double x) {
    return x * x;
}

double Vpt(double x, double n, double a) {
    return -n / square(cosh(a * x));
}

MatrixXd f(double x) {
    double v0 = Vpt(x, 45., 1.);
    double v1 = Vpt(x, 39. / 2, 1. / 2);
    MatrixXd r(2, 2);
    r << v0 + v1, v0 - v1, v0 - v1, v0 + v1;
    return r;
}
/*
TEST_CASE("Test propagatePsi", "[matscs][propagatePsi]") {
    Matscs<double> scs(&f, 2, 0, 20, 64);
    double Es[] = {-64, -36, -30.25, -20.25, -16};
    double mid = 20. / 64 * 2;
    for (double E : Es) {
        MatrixXd left = scs.propagatePsi(E, MatrixXd::Zero(2, 2), 0, mid);
        MatrixXd right = scs.propagatePsi(E, MatrixXd::Zero(2, 2), 20, mid);
        REQUIRE(Approx(0).margin(1e-3) == (left - right).determinant());
    }
}

*/