#include "../catch.hpp"
#include "../../src/matslise3d.h"
#include <iostream>

using namespace matslise;
using namespace std;

TEST_CASE("2D: coulomb", "[matslise3d][coulomb]") {
    // https://arxiv.org/pdf/0904.0939.pdf

    const double a = 0.05;
    const double b = 12.3;
    int steps = 16;
#pragma omp parallel for
    for (int i = 0; i < steps; ++i) {
        const double z = -b + 2 * b * (i + .5) / steps;
        Matslise2D<double> matslise([a, z](double x, double y) -> double {
            double r = hypot(hypot(x, y), z);
            if (r < a)
                return -2 / a;
            return -2 / r;
        }, {{-b, b}, -b, b}, Options2<double>().tolerance(1e-6).N(12).nested(Options1<double>().tolerance(1e-7)));


        double E;
        Eigen::Index j, m;
        Eigen::Index runningIndex = 0;
        for (auto &iEm : matslise.eigenvaluesByIndex(0, 10)) {
            tie(j, E, m) = iEm;
            REQUIRE(!isnan(E));
            REQUIRE(m > 0);
            REQUIRE(j == runningIndex);
            runningIndex += m;
        }
    }
}
