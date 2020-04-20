#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "../../src/util/constants.h"
#include "test_problem.h"

using namespace matslise;
using namespace std;
using namespace Eigen;


TEST_CASE("coffey_evans beta=20", "[matslise][coffey_evans]") {
    const double B = 20;
    MatsliseHalf<double> ms([B](double x) -> double {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, constants<double>::PI / 2, 1e-8);

    testProblem(ms, Y<>::Dirichlet(), Y<>::Dirichlet(), {
            -8.639369922249e-14, 77.91619567714, 151.4627783465, 151.4632236577, 151.4636689884, 220.1542298353,
            283.0948146954, 283.2507437431, 283.4087354034, 339.3706656525, 380.0949155509, 385.6447796090,
            394.1303198988, 426.5246237841, 452.6311747507, 477.7105126091, 507.5356903666, 540.6338227685,
            575.8375904214, 613.2813295704, 652.9904570847, 694.8904368843, 738.9381495539, 785.1084322837,
            833.3807330072, 883.7384298621, 936.1683286674, 990.6597832611, 1047.204086284, 1105.794050195,
    }, 1e-6);
}


TEST_CASE("coffey_evans beta=30", "[matslise][coffey_evans]") {
    const double B = 30;
    MatsliseHalf<double> ms([B](double x) -> double {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, constants<double>::PI / 2, 1e-14);

    testProblem(ms, Y<>::Dirichlet(), Y<>::Dirichlet(), {
            -2.629219539841e-14, 117.9463076621, 231.6649292371, 231.6649293130, 231.6649293888, 340.8882998096,
            445.2830895824, 445.2831723067, 445.2832550313, 544.4183851494, 637.6822498740, 637.7043623417,
            637.7265023185, 724.2576813500, 800.8313121817, 802.4787986926, 804.2793051689, 868.9602228708,
            909.4810465074, 925.9773098347, 951.8788067966, 992.5372967308, 1032.302482638, 1073.450031436,
            1118.144724204, 1165.668024705, 1215.559092002, 1267.808726318, 1322.384891028, 1379.228026450
    }, 1e-12, 1e-6);
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