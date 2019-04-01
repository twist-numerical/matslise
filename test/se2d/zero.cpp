#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include <matslise/se2d.h>


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Eigenvalues V=0", "[se2d][eigenfunctions][zero]") {
    SEnD<2> p2(
            [](double x, double y) -> double {
                return 0;
            },
            {{-M_PI_2, M_PI_2}, -M_PI_2, M_PI_2},
            Options<2>().sectorCount(13).nested(Options<1>().sectorCount(13)));
    double eigenvalues[] = {2, 5, 8, 10, 13, 17, 18};

    for (const double &E : eigenvalues) {
        const pair<double, double> &error = p2.calculateError(E);
        REQUIRE(abs(error.first) < 1e-9);
    }
}
