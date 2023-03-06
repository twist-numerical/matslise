#include <catch2/catch_test_macros.hpp>
#include "../../matslise/util/polynomial.h"
#include "fEquals.h"

using namespace matslise;


TEST_CASE("Test Polynomial", "[util][polynomial]") {

    Polynomial<int, 3> a({1, 2, 3, 4});

    REQUIRE(a.derivative<0>() == a);
    REQUIRE(a.derivative() == Polynomial<int, 2>{{2, 6, 12}});
    REQUIRE(a.derivative<2>() == Polynomial<int, 1>{{6, 24}});

    REQUIRE(a.integral() == Polynomial<int, 4>{{0, 1, 1, 1, 1}});

    REQUIRE(a * 2 == Polynomial<int, 3>{2, 4, 6, 8});
    REQUIRE(3 * a == Polynomial<int, 3>{3, 6, 9, 12});
    REQUIRE((3 * a) / 3 == a);

    REQUIRE(a - 2 == Polynomial<int, 3>{-1, 2, 3, 4});
    REQUIRE(3 + a == Polynomial<int, 3>{4, 2, 3, 4});

    REQUIRE((3 * (a - 9)) / 2 == ((9 - a) * -3) / 2);


    Polynomial<int, 2> b({5, 4, 3});


    REQUIRE(a * b == Polynomial<int, 5>{{5, 14, 26, 38, 25, 12}});


    Polynomial<double, 5> c{{5, 14, 26, 38, 25, 12}};

    fEquals<double>([&c](double x) { return c.derivative<0>(x); }, c, -10, 10);
    fEquals<double>([&c](double x) { return c.derivative<1>(x); }, c.derivative(), -10, 10);
    fEquals<double>([&c](double x) { return c.derivative<2>(x); }, c.template derivative<2>(), -10, 10);
}
