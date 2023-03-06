#include <cmath>
#include <vector>
#include <tuple>
#include "../test.h"
#include "../../matslise/matslise.h"
#include "../../matslise/util/constants.h"
#include "test_problem.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Marletta", "[matslise][marletta]") {
    Matslise<double> marletta([](double x) {
        return 3 * (x - 31) / (4 * (x + 1) * (x + 4) * (x + 4));
    }, 0, 12, 1e-9);
    Y<double> left({-8, 5}, {0, 0});
    Y<double> right({0, 1}, {0, 0});

    testProblem<double>(marletta, left, right, {
            -1.185214104750602, 0.044824295667483, 0.272327486048728, 0.657583587648663, 1.190706288916843,
            1.866791009935514, 2.683289140760672, 3.638839063660400, 4.732685922750347, 5.964393862874594,
            7.333701003760564, 8.840444083911645, 10.484517866821465, 12.265852456981172, 14.184400172534373,
            16.240127697625582, 18.433011249289144, 20.763033519061526, 23.230181689275661, 25.834446117162997,
            28.575819443870444
    }, 1e-6);
}

TEST_CASE("Marletta mirrored", "[matslise][marletta]") {
    Matslise<double> marletta([](double x) {
        x = 12 - x;
        return 3 * (x - 31) / (4 * (x + 1) * (x + 4) * (x + 4));
    }, 0, 12, 1e-9);
    Y<double> left({0, 1}, {0, 0});
    Y<double> right({8, 5}, {0, 0});

    testProblem<double>(marletta, left, right, {
            -1.185214104750602, 0.044824295667483, 0.272327486048728, 0.657583587648663, 1.190706288916843,
            1.866791009935514, 2.683289140760672, 3.638839063660400, 4.732685922750347, 5.964393862874594,
            7.333701003760564, 8.840444083911645, 10.484517866821465, 12.265852456981172, 14.184400172534373,
            16.240127697625582, 18.433011249289144, 20.763033519061526, 23.230181689275661, 25.834446117162997,
            28.575819443870444
    }, 1e-6);
}

TEST_CASE("Marletta mirrored boundary", "[matslise][marletta]") {
    Matslise<double> marletta([](double x) {
        return 3 * (x - 31) / (4 * (x + 1) * (x + 4) * (x + 4));
    }, 0, 12, 1e-9);
    Y<double> left({0, 1}, {0, 0});
    Y<double> right({-8, 5}, {0, 0});

    testProblem<double>(marletta, left, right, {
            -0.000630667376746, 0.148302684201381, 0.423712575204632, 0.833426710999476, 1.379612134380914,
            2.062996474420172, 2.883760193690856, 3.841907319455573, 4.937394289517807, 6.170171563439838,
            7.540195622512034, 9.047431228797695, 10.691850716010341, 12.473432593364979, 14.392160190135325,
            16.448020532712615, 18.641003468254635, 20.971100997740745, 23.438306773297711, 26.042615720079514,
            28.784023751215766
    }, 1e-6);
}
