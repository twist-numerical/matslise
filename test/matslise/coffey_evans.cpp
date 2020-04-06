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
                         const vector<pair<int, Scalar>> &eigenvalues, const Scalar &tolerance = Scalar(1e-7)) {
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
            REQUIRE(Approx(fx.y[0] * scale).margin(tolerance) == ys[j].y[0]);
            REQUIRE(Approx(fx.y[1] * scale).margin(tolerance) == ys[j].y[1]);
        }
    }
}

TEST_CASE("coffey_evans", "[matslise][coffey_evans]") {
    const double B = 20;
    Matslise<double> ms([B](double x) -> double {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, 0, constants<double>::PI / 2, 1e-9, sector_builder::uniform<Matslise<>>(31));

    Y<double> y0({0, 1}, {0, 0});
    Y<double> y1({1, 0}, {0, 0});
    vector<pair<int, double>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}

TEST_CASE("high potential", "[matslise][high]") {
    Matslise<double> ms([](double x) -> double {
        return (1 - cos(2 * constants<double>::PI * x)) / 2 * 1000;
    }, 0, 1, 1e-8);

    Y<double> y0({1, 0}, {0, 0});
    Y<double> y1({0, -1}, {0, 0});
    test_eigenfunctions(ms, y0, y1, ms.eigenvaluesByIndex(0, 20, y0, y1));
}

TEST_CASE("high potential (auto)", "[matslise][high]") {
    Matslise<double> ms([](double x) -> double {
        return (1 - cos(2 * constants<double>::PI * x)) / 2 * 1000;
    }, 0, 1, 1e-6);

    Y<double> y0({1, 0}, {0, 0});
    Y<double> y1({0, -1}, {0, 0});
    test_eigenfunctions(ms, y0, y1, ms.eigenvaluesByIndex(0, 20, y0, y1));
}

#ifdef MATSLISE_long_double

TEST_CASE("coffey_evans (long)", "[matslise][coffey_evans][long]") {
    const long double B = 20;
    Matslise<long double> ms([B](long double x) -> long double {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, 0, constants<long double>::PI / 2, 1e-9, sector_builder::uniform<Matslise<long double>>(31));

    Y<long double> y0({0, 1}, {0, 0});
    Y<long double> y1({1, 0}, {0, 0});
    vector<pair<int, long double>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}

TEST_CASE("high potential (long)", "[matslise][high][long]") {
    Matslise<long double> ms([](long double x) -> long double {
        return (1 - cos(2 * constants<double>::PI * x)) / 2 * 1000;
    }, 0, 1, 1e-9, sector_builder::uniform<Matslise<long double>>(31));

    Y<long double> y0({1, 0}, {0, 0});
    Y<long double> y1({0, -1}, {0, 0});
    vector<pair<int, long double>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}

TEST_CASE("high potential (auto) (long)", "[matslise][high][long]") {
    Matslise<long double> ms([](double x) -> double {
        return (1 - cos(2 * constants<double>::PI * x)) / 2 * 1000;
    }, 0, 1, 1e-11l);

    Y<long double> y0({1, 0}, {0, 0});
    Y<long double> y1({0, -1}, {0, 0});
    vector<pair<int, long double>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}

#endif

#ifdef MATSLISE_float128

#include <boost/multiprecision/float128.hpp>

using boost::multiprecision::float128;

TEST_CASE("coffey_evans (uniform)(float128)", "[matslise][coffey_evans][float128]") {
    const float128 B = 20;
    Matslise<float128> ms([B](float128 x) -> float128 {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, 0, constants<double>::PI, 1e-9q, sector_builder::uniform<Matslise<float128>>(31));

    Y<float128> y0({0, 1}, {0, 0});
    Y<float128> y1({1, 0}, {0, 0});
    vector<pair<int, float128>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}

TEST_CASE("high potential (float128)", "[matslise][high][float128]") {
    Matslise<float128> ms([](float128 x) -> float128 {
        return (1 - cos(2 * constants<double>::PI * x)) / 2 * 1000;
    }, 0, 1, 1e-16q, sector_builder::uniform<Matslise<float128>>(31));

    Y<float128> y0({1, 0}, {0, 0});
    Y<float128> y1({0, -1}, {0, 0});
    vector<pair<int, float128>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}

TEST_CASE("high potential (auto) (float128)", "[matslise][high][float128]") {
    Matslise<float128> ms([](float128 x) -> float128 {
        return (1 - cos(2 * constants<double>::PI * x)) / 2 * 1000;
    }, 0, 1, 1e-20q);

    Y<float128> y0({1, 0}, {0, 0});
    Y<float128> y1({0, -1}, {0, 0});
    vector<pair<int, float128>> eigenvalues = ms.eigenvaluesByIndex(0, 20, y0, y1);
    test_eigenfunctions(ms, y0, y1, eigenvalues);
}

#endif