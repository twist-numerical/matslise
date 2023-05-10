#include <cmath>
#include <iostream>
#include "../test.h"
#include "../../matslise/liouville.h"
#include "../../matslise/matslise.h"
#include "../../matslise/util/constants.h"


using namespace matslise;
using namespace Eigen;

inline double square(double x) {
    return x * x;
}

TEST_CASE("Klotter with liouville transformation", "[liouville]") {
    double a = 8. / 7;
    double b = 8;

    LiouvilleTransformation<double> transformation{
            {a, b},
            [](double) { return 1; },
            [](double x) { return 3 / (4 * x * x); },
            [](double x) { return square(8 * constants<double>::PI / (3 * x * x * x)); }
    };

    // std::cout << "Pieces: " << transformation.pieces.size() << std::endl;

    Matslise<double> m{[&](double x) { return transformation.V(x); }, transformation.xDomain().min,
                       transformation.xDomain().max};

    auto eigenvalues = m.eigenvaluesByIndex(0, 20, Y<double>::Dirichlet());
    Eigen::Index i = 0;
    for (auto &iE: eigenvalues) {
        double exact = square(iE.first + 1);
        REQUIRE(i == iE.first);
        REQUIRE_THAT(iE.second, WithinRel(exact, 1e-6));
        ++i;
        // std::cout << "E" << iE.first << ": " << iE.second << ", error: " << (iE.second - exact) << std::endl;
    }

}

TEST_CASE("Klotter as Sturm-Liouville problem", "[liouville]") {
    double a = 8. / 7;
    double b = 8;

    SturmLiouville<double> sl{
            [](double) { return 1; },
            [](double x) { return 3 / (4 * x * x); },
            [](double x) { return square(8 * constants<double>::PI / (3 * x * x * x)); },
            {a, b},
            1e-8
    };

    auto eigenvalues = sl.eigenvaluesByIndex(0, 20, Y<double>::Dirichlet(), Y<double>::Dirichlet());
    Eigen::Index i = 0;
    for (auto &iE: eigenvalues) {
        double exact = square(iE.first + 1);
        REQUIRE(i == iE.first);
        REQUIRE_THAT(iE.second, WithinRel(exact, 1e-6));
        ++i;
        // std::cout << "E" << iE.first << ": " << iE.second << ", error: " << (iE.second - exact) << std::endl;
    }

}
