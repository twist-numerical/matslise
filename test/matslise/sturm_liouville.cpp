#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../matslise/matslise.h"
#include "../../matslise/liouville.h"
#include "../../matslise/util/constants.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

// https://doi.org/10.17512/jamcm.2016.2.14

TEST_CASE("Example 1: p=1, q=0, w=(1+x)^-2", "[matslise][sturm-liouville]") {
    SturmLiouville<double> sl([](double) -> double { return 1; },
                              [](double) -> double { return 0; },
                              [](double x) { return 1 / (1 + x) / (1 + x); },
                              {0., 1.}, 1e-8);

    Y<> dirichlet = Y<>::Dirichlet();

    int i = 0;
    double f = constants<double>::PI / std::log(2);
    f *= f;
    for (auto &iE: sl.eigenvaluesByIndex(0, 20, dirichlet, dirichlet)) {
        REQUIRE(iE.first == i);
        double exact = 0.25 + f * (i + 1) * (i + 1);

        REQUIRE(Approx(exact).epsilon(1e-6) == iE.second);

        ++i;
    }
    REQUIRE(i == 20);
}

TEST_CASE("Example 2: p=1+x², q=x²-2, w=exp(x)", "[matslise][sturm-liouville]") {
    // The values in the article are not accurate. This true values are from Matslise 2.0.

    SturmLiouville<double> sl([](double x) { return (1 + x) * (1 + x); },
                              [](double x) { return x * x - 2; },
                              [](double x) { return std::exp(x); },
                              {0., 1.}, 1e-8);

    Y<> dirichlet = Y<>::Dirichlet();
    Y<> neumann = Y<>::Neumann();

    std::vector<double> exact
            {{
                     1.1704927599022, 26.8633676464100, 78.5490452641401, 156.0801562248233, 259.4554762109302,
                     388.6747906007730, 543.7380371293273, 724.6451922822132, 931.3962455677519, 1163.9911917274574,
                     1422.4300278884762
             }};
    int i = 0;
    for (auto &iE: sl.eigenvaluesByIndex(0, (int) exact.size(), dirichlet, neumann)) {
        REQUIRE(iE.first == i);

        REQUIRE(Approx(exact[i]).epsilon(1e-6) == iE.second);

        ++i;
    }
    REQUIRE(i == exact.size());
}

TEST_CASE("Example 3: p=2+sin(2πx), q=-10, w=1+sqrt(x)", "[matslise][sturm-liouville]") {
    // The values in the article are not accurate. This true values are from Matslise 2.0.
    // The stated boundary conditions are not the same as used in table 5
    // β_1 = 10, β_2 = 1

    SturmLiouville<double> sl([](double x) { return 2 + sin(2 * constants<double>::PI * x); },
                              [](double) -> double { return -10; },
                              [](double x) { return 1 + std::sqrt(x); },
                              {0., 1.}, 1e-8);

    Y<> dirichlet = Y<>::Dirichlet();
    Y<> right;
    right.y() << 2, -10;
    // a z + b p z' = 0, a = 10, b = 1
    // z_r = b p = 2
    // z_r' = -a = -10

    std::vector<double> exact
            {{
                     2.8978512448497, 24.0804524555819, 66.9259933904942, 130.2115116461903, 214.2598167805116,
                     319.2224520213227, 445.1429828280686, 592.0396118279557, 759.9206127332327, 948.7900946287176,
                     1158.6502115312053
             }};
    int i = 0;
    for (auto &iE: sl.eigenvaluesByIndex(0, (int) exact.size(), dirichlet, right)) {
        REQUIRE(iE.first == i);

        REQUIRE(Approx(exact[i]).epsilon(1e-6) == iE.second);

        ++i;
    }
    REQUIRE(i == exact.size());
}




