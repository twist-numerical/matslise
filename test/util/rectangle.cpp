#include "../catch.hpp"
#include "../../src/util/rectangle.h"

using namespace matslise;

TEST_CASE("rectangle", "[util][rectangle]") {
    Rectangle<double, 3> r{-3., 10., -11., -1., 0., 30.};
    REQUIRE(r.slice(0, 1, 2) == r.template slice<0, 1, 2>());
    REQUIRE(r.slice(2, 1, 0) == r.template slice<2, 1, 0>());
    REQUIRE(r.slice(0, 1, 2) != r.template slice<2, 1, 0>());
    REQUIRE(r.slice(0, 1, 0) == r.template slice<0, 1, 0>());

    REQUIRE(r.min(0) == r.min<0>());
    REQUIRE(r.min(0) == -3);
    REQUIRE(r.max(0) == r.max<0>());
    REQUIRE(r.max(0) == 10);

    REQUIRE(r.min(1) == r.min<1>());
    REQUIRE(r.min(1) == -11);
    REQUIRE(r.max(1) == r.max<1>());
    REQUIRE(r.max(1) == -1);

    REQUIRE(r.min(2) == r.min<2>());
    REQUIRE(r.min(2) == 0);
    REQUIRE(r.max(2) == r.max<2>());
    REQUIRE(r.max(2) == 30);
}
