#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
void test_eigenfunctions(const Matslise<Scalar> &ms, const Y<Scalar> &y0, const Y<Scalar> &y1,
                         const vector<pair<int, Scalar>> &eigenvalues) {
    Array<Scalar, Dynamic, 1> xs(100);
    int scalePoint = 39;
    for (int j = 0; j < xs.size(); ++j)
        xs[j] = ms.xmin + static_cast<Scalar>(j + .5) / static_cast<Scalar>(xs.size()) * (ms.xmax - ms.xmin);

    for (pair<int, Scalar> p : eigenvalues) {
        int i;
        Scalar E;
        tie(i, E) = p;
        std::function<Y<Scalar>(Scalar)> f = ms.eigenfunctionCalculator(E, y0, y1);
        Array<Y<Scalar>, Dynamic, 1> ys = ms.eigenfunction(E, y0, y1, xs);

        Scalar scale = ys[scalePoint].y[0] / f(xs[scalePoint]).y[0];
        for (int j = 0; j < xs.size(); ++j) {
            Y<Scalar> fx = f(xs[j]);
            REQUIRE(Approx(fx.y[0] * scale).margin(1e-7) == ys[j].y[0]);
            REQUIRE(Approx(fx.y[1] * scale).margin(1e-7) == ys[j].y[1]);
        }
    }
}

TEST_CASE("ixaru 1d", "[matslise][coffey_evans]") {
    Matslise<double> ms([](double x) -> double {
        const double y = 0;
        return (1 + x * x) * (1 + y * y);
    }, -5.5, 5.5, Matslise<double>::AUTO(1e-5));

    Y<double> y0({0, 1}, {0, 0});
    Y<double> y1({0, 1}, {0, 0});
    vector<pair<int, double>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    REQUIRE(eigenvalues.size() == 20);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}