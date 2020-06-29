#include "../catch.hpp"
#include "../../src/matslise3d.h"
#include <iostream>

using namespace matslise;
using namespace std;

TEST_CASE("3D: coulomb", "[matslise3d][coulomb]") {
    // https://arxiv.org/pdf/0904.0939.pdf

    const double a = 0.05;
    Matslise3D<double> matslise([a](double x, double y, double z) -> double {
        double r = hypot(hypot(x, y), z);
        if (r < a)
            return -1 / a;
        return -1 / r;
    }, {{{-12.3, 12.3}, -12.3, 12.3}, -12.3, 12.3}, sector_builder::uniform<Matslise3D<double>>(18), 1e-7);

    cout << matslise.eigenvalue(-0.48) << endl;
    cout << matslise.eigenvalue(-0.50) << endl;
    cout << matslise.eigenvalue(-0.52) << endl;
    cout << endl;
    cout << matslise.eigenvalue(-0.11) << endl;
    cout << matslise.eigenvalue(-0.122) << endl;
    cout << matslise.eigenvalue(-0.13) << endl;
}
