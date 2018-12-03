#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include <matslise/matslise.h>


using namespace matslise;
using namespace std;
using namespace Eigen;

double mathieu(double x) {
    return 2 * cos(2 * x);
}

TEST_CASE("Solving the mathieu problem (first 10)", "[matslise][mathieu]") {
    Matslise ms(&mathieu, 0, M_PI, 8);

    vector<double> correct = {-0.11024881635796, 3.91702477214389, 9.04773925867679, 16.03297008079835,
                              25.02084082368434, 36.01428991115492, 49.01041825048373, 64.00793719066102,
                              81.00625032615399, 100.00505067428990, 121.00416676119610};

    vector<pair<int, double>> *eigenvalues = ms.computeEigenvaluesByIndex(0, (int) correct.size(),
                                                                                   Y<double>({0, 1}),
                                                                                   Y<double>({0, 1}));
    for (unsigned int i = 0; i < correct.size(); ++i) {
        REQUIRE(i == eigenvalues->at(i).first);
        REQUIRE(Approx(correct[i]).margin(1e-12) == eigenvalues->at(i).second);
    }
    delete eigenvalues;
}

TEST_CASE("Solving the mathieu problem (skip 100)", "[matslise][mathieu]") {
    Matslise ms(&mathieu, 0, M_PI, 8);

    vector<double> correct = {10201.000049019607, 10404.000048063059, 10609.000047134254, 10816.000046232084,
                              11025.000045355584, 11236.000044503782, 11449.000043675751, 11664.000042870637,
                              11881.000042087542, 12100.000041325728, 12321.000040584415};
    unsigned int offset = 100;
    vector<pair<int, double>> *eigenvalues = ms.computeEigenvaluesByIndex(
            offset, offset + (unsigned int) correct.size(), Y<double>({0, 1}), Y<double>({0, 1}));

    REQUIRE(correct.size() == eigenvalues->size());
    for (unsigned int i = 0; i < correct.size(); ++i) {
        REQUIRE(offset + i == eigenvalues->at(i).first);
        REQUIRE(Approx(correct[i]).margin(1e-12) == eigenvalues->at(i).second);
    }
    delete eigenvalues;
}

TEST_CASE("Mathieu problem eigenfunctions", "[mathieu][matslise][eigenfunctions]") {
    ArrayXd x(30);
    x << 0.000000000000000, 0.108330781158269, 0.216661562316537, 0.324992343474806, 0.433323124633075,
                        0.541653905791344, 0.649984686949612, 0.758315468107881, 0.866646249266150, 0.974977030424419,
                        1.083307811582687, 1.191638592740956, 1.299969373899225, 1.408300155057494, 1.516630936215762,
                        1.624961717374031, 1.733292498532300, 1.841623279690568, 1.949954060848837, 2.058284842007106,
                        2.166615623165375, 2.274946404323643, 2.383277185481912, 2.491607966640180, 2.599938747798450,
                        2.708269528956718, 2.816600310114987, 2.924931091273256, 3.033261872431524, 3.141592653589793;

    vector<double> y0 = {-0.00000000000000, 0.05957648373348, 0.12058229947745, 0.18428570118601, 0.25163008589212,
                         0.32306954040196, 0.39841414637477, 0.47670540198060, 0.55615034750173, 0.63414569796321,
                         0.70741709118712, 0.77228194717948, 0.82501918024869, 0.86230069024582, 0.88161654418474,
                         0.88161654418474, 0.86230069024582, 0.82501918024869, 0.77228194717948, 0.70741709118712,
                         0.63414569796321, 0.55615034750173, 0.47670540198060, 0.39841414637477, 0.32306954040196,
                         0.25163008589212, 0.1842857016118601, 0.12058229947745, 0.05957648373348, 0.00000000000000};

    vector<double> dy0 = {0.54770139190616, 0.55442179455968, 0.57384238641313, 0.60372505966206, 0.64030804382969,
                          0.67833289106582, 0.71121342385303, 0.73143606944606, 0.73124805111186, 0.70362163948128,
                          0.64338934196303, 0.54835108763979, 0.42009484461026, 0.26427976960890, 0.09022152372952,
                          -0.09022152372952, -0.26427976960890, -0.42009484461026, -0.54835108763979, -0.64338934196303,
                          -0.70362163948128, -0.73124805111186, -0.73143606944606, -0.71121342385303, -0.67833289106582,
                          -0.64030804382969, -0.60372505966206, -0.57384238641313, -0.55442179455968,
                          -0.54770139190616};

    vector<double> y3 = {-0.00000000000000, 0.32417284029898, 0.59543237414050, 0.76858796800638, 0.81320927029333,
                         0.71916231964893, 0.49966578556769, 0.19086652240175, -0.15277306606069, -0.46727765184401,
                         -0.69110153718396, -0.77802995899298, -0.70806542774819, -0.49329149619403, -0.17665033540715,
                         0.17665033540715, 0.49329149619403, 0.70806542774819, 0.77802995899298, 0.69110153718396,
                         0.46727765184401, 0.15277306606069, -0.19086652240175, -0.49966578556769, -0.71916231964893,
                         -0.81320927029333, -0.76858796800638, -0.59543237414050, -0.32417284029898, -0.00000000000000};


    vector<double> dy3 = {3.07626133037988, 2.82599938325101, 2.11136070876314, 1.03774249508994, -0.23002968436917,
                          -1.48587759669341, -2.51014152163739, -3.10573501735402, -3.13821749458431, -2.57141674375545,
                          -1.48796775591086, -0.08518552873415, 1.35809781500210, 2.53962297948860, 3.20361729625403,
                          3.20361729625403, 2.53962297948859, 1.35809781500210, -0.08518552873415, -1.48796775591086,
                          -2.57141674375544, -3.13821749458430, -3.10573501735402, -2.51014152163740, -1.48587759669341,
                          -0.23002968436916, 1.03774249508994, 2.11136070876314, 2.82599938325101, 3.07626133037988};

    Matslise ms(&mathieu, 0, M_PI, 8);
    Y<double> ystart({0,1});

    {
        auto *eigenvalues = ms.computeEigenvaluesByIndex(0, 1, ystart, ystart);
        REQUIRE(eigenvalues->size() == 1);
        REQUIRE(0 == eigenvalues->at(0).first);
        double e = eigenvalues->at(0).second;
        delete eigenvalues;

        REQUIRE(Approx(-0.11024881635796).margin(1e-12) == e);
        Array<Y<double>, Dynamic, 1> result = ms.computeEigenfunction(e, ystart, ystart, x);
        REQUIRE(result.size() == y0.size());
        for (long i = result.size() - 1; i >= 0; --i)
            result[i] *= dy0[0] / result[0].y[1];

        for (int i = 0; i < result.size(); ++i) {
            REQUIRE(Approx(y0[i]).margin(1e-12) == result[i].y[0]);
            REQUIRE(Approx(dy0[i]).margin(1e-12) == result[i].y[1]);
        }
    }

    {
        auto *eigenvalues = ms.computeEigenvaluesByIndex(3, 4, ystart, ystart);
        REQUIRE(eigenvalues->size() == 1);
        double e = eigenvalues->at(0).second;
        delete eigenvalues;

        REQUIRE(Approx(16.03297008140580).margin(1e-12) == e);
        Array<Y<double>, Dynamic, 1> result = ms.computeEigenfunction(e, ystart, ystart, x);
        REQUIRE(result.size() == y0.size());
        for (long i = result.size() - 1; i >= 0; --i)
            result[i] *= dy3[0] / result[0].y[1];

        for (int i = 0; i < result.size(); ++i) {
            REQUIRE(Approx(y3[i]).margin(1e-12) == result[i].y[0]);
            REQUIRE(Approx(dy3[i]).margin(1e-12) == result[i].y[1]);
        }
    }
}