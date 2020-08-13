#include "../catch.hpp"
#include "../../src/matslise3d.h"
#include <iostream>

using namespace matslise;
using namespace std;

TEST_CASE("3D: coulomb", "[matslise3d][coulomb]") {
    // https://arxiv.org/pdf/0904.0939.pdf

    Matslise3D<>::Config config;
    config.tolerance = 1e-7;
    config.xyBasisSize = 14;
    config.zStepsPerSector = 3;
    config.xSymmetric = true;
    config.ySymmetric = true;
    config.zSectorBuilder = sector_builder::uniform<Matslise3D<>>(16);
    config.ySectorBuilder = sector_builder::automatic<Matslise2D<>>(1e-7);

    const double a = 0.01;
    Matslise3D<double> matslise([a](double x, double y, double z) -> double {
        double r = hypot(hypot(x, y), z);
        if (r < a)
            return -2 / a;
        return -2 / r;
    }, {{{-12.3, 12.3}, -12.3, 12.3}, -12.3, 12.3}, config);

    for (double Eguess : (double[]) {-1, -.25, -.1111}) {
        double E = matslise.eigenvalue(Eguess);
        cout << E << " (±" << matslise.eigenvalueError(E) << ")" << endl;
    }
}
