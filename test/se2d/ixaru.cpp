#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Eigenfunctions ixaru", "[se2d][eigenfunctions][ixaru]") {
    SEnD<2> p2(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {{-5.5, 5.5}, -5.5, 5.5},
            Options<2>().sectorCount(23).stepsPerSector(4).N(10).nested(Options<1>().sectorCount(26)));
    pair<double, int> eigenvalues[] = {
            {3.1959181,  1},
            {5.5267439,  2},
            {7.5578033,  1},
            {8.0312723,  1},
            {8.4445814,  1},
            {9.9280611,  2},
            {11.3118171, 2},
            {12.1032536, 1},
            {12.2011790, 1},
            {13.3323313, 1}
    };


    ArrayXd x(3);
    x << -1, 0, 1;
    double E;
    int multiplicity;
    for (int i = 0; i < 9; ++i) {
        double error = get<0>(p2.calculateError((get<0>(eigenvalues[i]) + get<0>(eigenvalues[i + 1])) / 2));
        CHECK(abs(error) > 1e-3);
    }
    vector<double> eigenvalues_simple;
    for (auto &Emult  : eigenvalues) {
        tie(E, multiplicity) = Emult;
        eigenvalues_simple.push_back(E);
        const pair<double, double> &error = p2.calculateError(E);
        CHECK(abs(error.first) < 1e-3);

        if (E < 6) {
            double El = p2.findEigenvalue(E - 0.01);
            double Em = p2.findEigenvalue(E + 0.01);
            CHECK(Approx(El).margin(1e-7) == E);
            CHECK(Approx(Em).margin(1e-7) == E);
        }

        const vector<Array<double, -1, -1>> *f = p2.computeEigenfunction(E, {x, x});
        CHECK(f->size() == multiplicity);
        delete f;
    }
    checkOrthonormality(p2, eigenvalues_simple.begin(), eigenvalues_simple.end());
}

TEST_CASE("Eigenfunctions ixaru auto", "[se2d][eigenfunctions][ixaru][auto]") {
    SEnD<2> p2(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {{-5.5, 5.5}, -5.5, 5.5},
            Options<2>().tolerance(1e-5).stepsPerSector(3).N(10).nested(Options<1>().tolerance(1e-7)));
    pair<double, int> eigenvalues[] = {
            {3.1959181,  1},
            {5.5267439,  2},
            {7.5578033,  1},
            {8.0312723,  1},
            {8.4445814,  1},
            {9.9280611,  2},
            {11.3118171, 2},
            {12.1032536, 1},
            {12.2011790, 1},
            {13.3323313, 1}
    };


    ArrayXd x(3);
    x << -1, 0, 1;
    double E;
    int multiplicity;
    for (int i = 0; i < 9; ++i) {
        double error = get<0>(p2.calculateError((get<0>(eigenvalues[i]) + get<0>(eigenvalues[i + 1])) / 2));
        CHECK(abs(error) > 1e-3);
    }
    vector<double> eigenvalues_simple;
    for (auto &Emult  : eigenvalues) {
        tie(E, multiplicity) = Emult;
        eigenvalues_simple.push_back(E);
        const pair<double, double> &error = p2.calculateError(E);
        CHECK(abs(error.first) < 1e-3);

        if (E < 6) {
            double El = p2.findEigenvalue(E - 0.01);
            double Em = p2.findEigenvalue(E + 0.01);
            CHECK(Approx(El).margin(1e-7) == E);
            CHECK(Approx(Em).margin(1e-7) == E);
        }

        const vector<Array<double, -1, -1>> *f = p2.computeEigenfunction(E, {x, x});
        CHECK(f->size() == multiplicity);
        delete f;
    }
    checkOrthonormality(p2, eigenvalues_simple.begin(), eigenvalues_simple.end());
}
