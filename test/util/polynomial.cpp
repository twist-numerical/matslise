#include "../catch.hpp"
#include "../../matslise/util/polynomial.h"

using namespace matslise;


TEST_CASE("Test Polynomial", "[util][polynomial]") {

    Polynomial<int, 3> a({1, 2, 3, 4});

    REQUIRE(a.derivative<0>() == a);
    REQUIRE(a.derivative() == Polynomial<int, 2>{{2, 6, 12}});
    REQUIRE(a.derivative<2>() == Polynomial<int, 1>{{6, 24}});

    REQUIRE(a.integral() == Polynomial<int, 4>{{0, 1, 1, 1, 1}});


    Polynomial<int, 2> b({5, 4, 3});


    REQUIRE(a * b == Polynomial<int, 5>{{5, 14, 26, 38, 25, 12}});


    Polynomial<double, 5> c{{5, 14, 26, 38, 25, 12}};

    auto dc = c.derivative();
    auto ddc = c.derivative<2>();
    for (double x = -10; x <= 10; ++x) {
        REQUIRE(Approx(c.derivative<0>(x)) == c(x));
        REQUIRE(Approx(c.derivative<1>(x)) == dc(x));
        REQUIRE(Approx(c.derivative<2>(x)) == ddc(x));
    }

}
