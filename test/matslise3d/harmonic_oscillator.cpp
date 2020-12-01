#include "../catch.hpp"
#include "../../src/matslise3d.h"
#include <iostream>

using namespace matslise;
using namespace std;

/* IGNORING FOR NOW
TEST_CASE("3D: harmonic oscillator", "[matslise3d][harmonic][slow]") {
    // https://arxiv.org/pdf/0904.0939.pdf

    Matslise3D<>::Config config;
    config.tolerance = 1e-8;
    config.zStepsPerSector = 2;
    config.yStepsPerSector = 2;
    config.xBasisSize = 18;
    config.xyBasisSize = 22;
    config.zSectorBuilder = sector_builder::automatic<Matslise3D<>, true>(1e-6);

    const double a = 0.01;
    Matslise3D<double> matslise([](double x, double y, double z) -> double {
        return hypot(hypot(x, y), z);
    }, {-6., 6., -6., 6., -6., 6.}, config);

    for (auto &iEm : matslise.eigenvaluesByIndex(0, 10)) {
        cout << get<0>(iEm) << ", " << get<1>(iEm) << ", " << get<2>(iEm) << endl;
    }
}
*/
