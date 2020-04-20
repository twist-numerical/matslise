#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "test_problem.h"


using namespace matslise;
using namespace std;
using namespace Eigen;


TEST_CASE("ixaru 1d", "[matslise][coffey_evans]") {
    Matslise<double> ms([](double x) -> double {
        const double y = 0;
        return (1 + x * x) * (1 + y * y);
    }, -5.5, 5.5, 1e-5);

    Y<double> y0({0, 1}, {0, 0});
    Y<double> y1({0, 1}, {0, 0});
    vector<pair<int, double>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    REQUIRE(eigenvalues.size() == 20);
    testOrthogonality<double>(ms, y0, y1, eigenvalues, 1e-4);
}