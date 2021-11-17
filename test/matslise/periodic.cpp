#include "../catch.hpp"
#include "../../matslise/matslise.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Periodic zero potential", "[matslise][periodic]") {
    PeriodicMatslise<> matslise([](double) { return 0; }, 0, M_PI, 1e-8);
    auto error = matslise.matchingError(3);

}

TEST_CASE("Periodic mathieu", "[matslise][periodic]") {
    PeriodicMatslise<> matslise([](double x) { return 2 * cos(2 * x); }, 0, M_PI, 1e-8);

    for (auto &eigenpair : matslise.eigenpairsByIndex(0, 20)) {
        cout << "i: " << get<0>(eigenpair) << ", E: " << get<1>(eigenpair) << ", m: " << get<2>(eigenpair).size()
             << endl;
    }
}

