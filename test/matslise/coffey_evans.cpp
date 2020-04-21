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

TEST_CASE("coffey_evans beta=20 (long)", "[matslise][coffey_evans][long]") {
    const long double B = 20.l;
    MatsliseHalf<long double> ms([B](long double x) -> long double {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, constants<long double>::PI / 2, 1e-12l);

    testProblem<long double>(ms, Y<long double>::Dirichlet(), Y<long double>::Dirichlet(), {
            -0.00000000000002l, 77.91619567714397l, 151.46277834645664l, 151.46322365765860l, 151.46366898835166l,
            220.15422983525997l, 283.09481469540145l, 283.25074374311265l, 283.40873540342932l, 339.37066565252246l,
            380.09491555093172l, 385.64477960900808l, 394.13031989879846l, 426.52462378409643l, 452.63117475070641l,
            477.71051260907672l, 507.53569036662458l, 540.63382276850348l, 575.83759042140582l, 613.28132957039725l,
            652.99045708465667l
    }, 1e-10l, 1e-6l);
}

TEST_CASE("coffey_evans beta=30 (long)", "[matslise][coffey_evans][long]") {
    const long double B = 30;
    MatsliseHalf<long double> ms([B](long double x) -> long double {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, constants<long double>::PI / 2, 1e-16l);

    testProblem<long double>(ms, Y<long double>::Dirichlet(), Y<long double>::Dirichlet(), {
            0.0l, 117.94630766207001l, 231.66492923712826l, 231.66492931296196l, 231.66492938879651l,
            340.88829980961299l, 445.28308958243815l, 445.28317230667329l, 445.28325503133368l, 544.41838514936001l,
            637.68224987405006l, 637.70436234165697l, 637.72650231851446l, 724.25768135002033l, 800.83131218170922l,
            802.47879869262465l, 804.27930516886272l, 868.96022287077210l, 909.48104650741504l, 925.97730983473230l,
            951.87880679659270l
    }, 1e-12l, 1e-6l);
}

#endif

#ifdef MATSLISE_float128

#include <boost/multiprecision/float128.hpp>

using boost::multiprecision::float128;

TEST_CASE("coffey_evans beta=20 (float128)", "[matslise][coffey_evans][float128]") {
    const float128 B = 20;
    MatsliseHalf<float128> ms([B](float128 x) -> float128 {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, constants<float128>::PI / 2, 1e-14q);

    testProblem<float128>(ms, Y<float128>::Dirichlet(), Y<float128>::Dirichlet(), {
            -0.00000000000002q, 77.91619567714397q, 151.46277834645664q, 151.46322365765860q, 151.46366898835166q,
            220.15422983525997q, 283.09481469540145q, 283.25074374311265q, 283.40873540342932q, 339.37066565252246q,
            380.09491555093172q, 385.64477960900808q, 394.13031989879846q, 426.52462378409643q, 452.63117475070641q,
            477.71051260907672q, 507.53569036662458q, 540.63382276850348q, 575.83759042140582q, 613.28132957039725q,
            652.99045708465667q
    }, 1e-12q, 1e-6q);
}

TEST_CASE("coffey_evans beta=30 (float128)", "[matslise][coffey_evans][float128]") {
    const float128 B = 30;
    MatsliseHalf<float128> ms([B](float128 x) -> float128 {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, constants<float128>::PI / 2, 1e-16q);

    testProblem<float128>(ms, Y<float128>::Dirichlet(), Y<float128>::Dirichlet(), {
            0.0q, 117.94630766207001q, 231.66492923712826q, 231.66492931296196q, 231.66492938879651q,
            340.88829980961299q, 445.28308958243815q, 445.28317230667329q, 445.28325503133368q, 544.41838514936001q,
            637.68224987405006q, 637.70436234165697q, 637.72650231851446q, 724.25768135002033q, 800.83131218170922q,
            802.47879869262465q, 804.27930516886272q, 868.96022287077210q, 909.48104650741504q, 925.97730983473230q,
            951.87880679659270q
    }, 1e-12q, 1e-6q);
}

#endif