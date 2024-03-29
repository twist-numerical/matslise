#include <cmath>
#include <vector>
#include <tuple>
#include "../test.h"
#include "../../matslise/matscs.h"


using namespace matslise;
using namespace matslise::sector_builder;
using namespace std;
using namespace Eigen;

TEST_CASE("Test propagateTheta", "[matscs][propagateTheta]") {
    Matscs<double> scs([](double x) -> MatrixXd {
        MatrixXd r(4, 4);
        double dx = 1;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
                r(i, j) = cos(x) / max(i + 1, j + 1);
            dx *= x;
            r(i, i) += 1 / dx;
        }
        return r;
    }, 4, 0.1, 1, AutomaticSectorBuilder<Matscs<>>(1e-5));
    double E = 14.94180054416473;
    double theta = 0;
    REQUIRE_THAT(scs.propagate(E, Y<double, Dynamic>::Dirichlet(4), 0.1, 1, theta).first
                         .block().determinant(), WithinAbs<double>(0, 1e-3));
}
