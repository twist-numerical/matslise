#include "../catch.hpp"
#include "../../matslise/matslise.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Periodic zero potential", "[matslise][periodic]") {
    PeriodicMatslise<> matslise([](double) { return 0; }, 0, M_PI, 1e-8);
    auto error = matslise.matchingError(3);
    cout << get<0>(error) << endl;
    cout << get<2>(error).transpose() << endl;
}

