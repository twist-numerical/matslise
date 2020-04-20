#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "../../src/util/constants.h"
#include "test_problem.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Marletta", "[matslise][marletta]") {
    Matslise<double> marletta([](double x) {
        return 3 * (x - 31) / (4 * (x + 1) * (x + 4) * (x + 4));
    }, 0, 12, 1e-6);
    Y<double> left({-8, 5}, {0, 0});
    Y<double> right({0, 1}, {0, 0});

    testProblem<double>(marletta, left, right, {
            0, 0.04482429566748, 0.2723274860487, 0.6575835876487, 1.190706288917, 1.866791009936, 2.683289140761,
            3.638839063660, 4.732685922750, 5.964393862875, 7.333701003761
    }, 1e-6);
}
