#include <cmath>
#include <iostream>
#include "../catch.hpp"
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
            [](double x) { return 1; },
            [](double x) { return 3 / (4 * x * x); },
            [](double x) { return square(8 * constants<double>::PI / (3 * x * x * x)); }
    };

    std::cout << "Pieces: " << transformation.pieces.size() << std::endl;

    Matslise<double> m{[&](double x) { return transformation.V(x); }, transformation.xDomain().min(),
                       transformation.xDomain().max()};

    auto eigenvalues = m.eigenvaluesByIndex(0, 20, Y<double>::Dirichlet());
    for (auto &iE : eigenvalues) {
        double exact = square(iE.first + 1);
        std::cout << "E" << iE.first << ": " << iE.second << ", error: " << (iE.second - exact) << std::endl;
    }

}
