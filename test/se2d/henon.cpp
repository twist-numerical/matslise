#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/se2d.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Eigenfunctions henon", "[se2d][eigenfunctions][henon]") {
    auto *p2 = new SEnD<2>(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2 * sqrt(5)) * x * (y * y - x * x / 3);
            },
            {{-6, 6}, -6, 6},
            Options<2>().sectorCount(16).stepsPerSector(4).nested(Options<1>().sectorCount(16)));
    pair<double, int> eigenvalues[] = {
            {0.998594690530479, 1},
            {1.99007660445524, 2},
            {2.95624333869018, 1},
            {2.98532593386986, 2},
            {3.92596412795287, 2},
            {3.98241882458866, 1},
            {3.98575763690663, 1},
            {4.87014557482289, 1},
            {4.89864497284387, 2}};

    ArrayXd x(3);
    x << -1, 0, 1;
    double E;
    int multiplicity;
    for (auto &Em  : eigenvalues) {
        tie(E, multiplicity) = Em;
        E *= 2;

        const pair<double, double> &error = p2->calculateError(E);
        REQUIRE(abs(error.first) < 1e-3);

        const vector<Array<double, -1, -1>> *f = p2->computeEigenfunction(E, {x, x});
        REQUIRE(f->size() == multiplicity);
        delete f;
    }
    delete p2;
}
