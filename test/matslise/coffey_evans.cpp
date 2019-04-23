#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

void test_eigenfunctions(const Matslise &ms, const Y<> &y0, const Y<> &y1, vector<pair<int, double>> *eigenvalues) {


    ArrayXd xs(100);
    int scalePoint = 39;
    for (int j = 0; j < xs.size(); ++j)
        xs[j] = ms.xmin + (j + .5) / xs.size() * (ms.xmax - ms.xmin);

    for (pair<int, double> p : *eigenvalues) {
        int i;
        double E;
        tie(i, E) = p;
        std::function<Y<>(double)> f = ms.eigenfunctionCalculator(E, y0, y1);
        Array<Y<>, Dynamic, 1> ys = ms.computeEigenfunction(E, y0, y1, xs);

        double scale = ys[scalePoint].y[0] / f(xs[scalePoint]).y[0];
        for (int j = 0; j < xs.size(); ++j) {
            Y<> fx = f(xs[j]);
            REQUIRE(Approx(fx.y[0] * scale).margin(1e-7) == ys[j].y[0]);
            REQUIRE(Approx(fx.y[1] * scale).margin(1e-7) == ys[j].y[1]);
        }
    }
}

TEST_CASE("coffey_evans", "[matslise][coffey_evans]") {
    const double B = 20;
    Matslise ms([B](double x) -> double {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, 0, M_PI_2, 31);

    Y<> y0({0, 1}, {0, 0});
    Y<> y1({1, 0}, {0, 0});
    vector<pair<int, double>> *eigenvalues = ms.computeEigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
    delete eigenvalues;
}

TEST_CASE("high potential", "[matslise][high]") {
    Matslise ms([](double x) -> double {
        return (1 - cos(2 * M_PI * x)) / 2 * 1000;
    }, 0, 1, 31);

    Y<> y0({1, 0}, {0, 0});
    Y<> y1({0, -1}, {0, 0});
    vector<pair<int, double>> *eigenvalues = ms.computeEigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
    delete eigenvalues;
}

TEST_CASE("high potential (auto)", "[matslise][high]") {
    Matslise ms([](double x) -> double {
        return (1 - cos(2 * M_PI * x)) / 2 * 1000;
    }, 0, 1, Matslise::AUTO(1e-6));

    Y<> y0({1, 0}, {0, 0});
    Y<> y1({0, -1}, {0, 0});
    vector<pair<int, double>> *eigenvalues = ms.computeEigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
    delete eigenvalues;
}