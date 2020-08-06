#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

TEST_CASE("Eigenfunctions ixaru", "[matslise2d][eigenfunctions][ixaru]") {
    Matslise2D<> p2(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {{-5, 5}, -5, 5},
            Options2<>().tolerance(1e-7).nested(Options1<>().tolerance(1e-7)));
    //Options2<>().sectorCount(23).stepsPerSector(4).N(10).nested(Options1<>().sectorCount(26)));
    pair<double, unsigned int> eigenvalues[] = {
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

    const vector<tuple<Index, double, Index>> foundEigenvalues = p2.eigenvaluesByIndex(0, 10);
    int j = 0;
    for (int i = 0; i < 10; ++i) {
        double E = get<1>(foundEigenvalues[i]);
        CHECK(j == get<0>(foundEigenvalues[i]));
        CHECK(eigenvalues[i].second == get<2>(foundEigenvalues[i]));
        CHECK(Approx(E).margin(1e-7) == eigenvalues[i].first);
        CHECK(p2.estimateIndex(E - 1e-3) == j);
        j += eigenvalues[i].second;
        CHECK(p2.estimateIndex(E + 1e-3) == j);
    }

    double E;
    unsigned int multiplicity;
    for (int i = 0; i < 9; ++i) {
        double error = get<0>(p2.matchingError((get<0>(eigenvalues[i]) + get<0>(eigenvalues[i + 1])) / 2));
        CHECK(abs(error) > 1e-3);
    }
    vector<double> eigenvalues_simple;
    for (auto &Emult  : eigenvalues) {
        tie(E, multiplicity) = Emult;
        eigenvalues_simple.push_back(E);
        const pair<double, double> &error = p2.matchingError(E);
        CHECK(abs(error.first) < 1e-3);
        CHECK(p2.eigenvalueError(E) < 1e-5);

        if (E < 6) {
            double El = p2.eigenvalue(E - 0.01).first;
            double Em = p2.eigenvalue(E + 0.01).first;
            CHECK(Approx(El).margin(1e-7) == E);
            CHECK(Approx(Em).margin(1e-7) == E);
        }

        const vector<Eigenfunction2D<>> f = p2.eigenfunction(E);
        CHECK(f.size() == multiplicity);
    }
    checkOrthonormality(&p2, eigenvalues_simple.begin(), eigenvalues_simple.end());
}

TEST_CASE("Eigenfunctions ixaru halfrange", "[matslise2d][eigenfunctions][ixaru][halfrange]") {
    Matslise2DHalf<> p2(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {{-5.5, 5.5}, -5.5, 5.5},
            Options2<>().tolerance(1e-6).N(10).nested(Options1<>().sectorCount(26)));
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

    // TODO: Some issues with matchingError function of half dirichlet-problem
    /*
    const vector<double> foundEigenvalues = p2.eigenvaluesByIndex(0, 10);
    for (int i = 0; i < 10; ++i) {
        CHECK(Approx(foundEigenvalues[i]).margin(1e-7) == eigenvalues[i].first);
    }
     */

    double E;
    int multiplicity;

    vector<double> eigenvalues_simple;
    for (auto &Emult  : eigenvalues) {
        tie(E, multiplicity) = Emult;
        eigenvalues_simple.push_back(E);

        if (E < 6) {
            double El = p2.eigenvalue(E - 0.01).first;
            double Em = p2.eigenvalue(E + 0.01).first;
            CHECK(Approx(El).margin(1e-7) == E);
            CHECK(Approx(Em).margin(1e-7) == E);
            CHECK(p2.eigenvalueError(El) < 1e-5);
            CHECK(p2.eigenvalueError(Em) < 1e-5);
        }

        const vector<Eigenfunction2D<>> f = p2.eigenfunction(E);
        CHECK(static_cast<long>(f.size()) == multiplicity);
    }
    checkOrthonormality(&p2, eigenvalues_simple.begin(), eigenvalues_simple.end());
}

TEST_CASE("Eigenfunctions ixaru auto", "[matslise2d][eigenfunctions][ixaru][auto]") {
    Matslise2D<> p2(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {{-5.5, 5.5}, -5.5, 5.5},
            Options2<>().tolerance(1e-6).N(10).nested(Options1<>().tolerance(1e-7)));
    pair<double, unsigned int> eigenvalues[] = {
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

    double E;
    unsigned int multiplicity;
    for (int i = 0; i < 9; ++i) {
        double error = get<0>(p2.matchingError((get<0>(eigenvalues[i]) + get<0>(eigenvalues[i + 1])) / 2));
        CHECK(abs(error) > 1e-3);
    }
    vector<double> eigenvalues_simple;
    for (auto &Emult  : eigenvalues) {
        tie(E, multiplicity) = Emult;
        eigenvalues_simple.push_back(E);
        const pair<double, double> &error = p2.matchingError(E);
        CHECK(abs(error.first) < 1e-3);

        if (E < 6) {
            double El = p2.eigenvalue(E - 0.01).first;
            double Em = p2.eigenvalue(E + 0.01).first;
            CHECK(Approx(El).margin(1e-7) == E);
            CHECK(Approx(Em).margin(1e-7) == E);
            CHECK(p2.eigenvalueError(El) < 1e-5);
            CHECK(p2.eigenvalueError(Em) < 1e-5);
        }

        const vector<Eigenfunction2D<>> f = p2.eigenfunction(E);
        CHECK(f.size() == multiplicity);
    }
    checkOrthonormality(&p2, eigenvalues_simple.begin(), eigenvalues_simple.end());
}

TEST_CASE("Eigenfunctions ixaru auto high n", "[matslise2d][eigenfunctions][ixaru][auto][slow]") {
    Matslise2D<> p2(
            [](double x, double y) -> double {
                return (1 + x * x) * (1 + y * y);
            },
            {{-5.5, 5.5}, -5.5, 5.5},
            Options2<>().tolerance(1e-5).N(20).nested(Options1<>().tolerance(1e-7)));
    pair<double, unsigned int> eigenvalues[] = {
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

    double E;
    unsigned int multiplicity;
    for (int i = 0; i < 9; ++i) {
        double error = get<0>(p2.matchingError((get<0>(eigenvalues[i]) + get<0>(eigenvalues[i + 1])) / 2));
        CHECK(abs(error) > 1e-3);
    }
    vector<double> eigenvalues_simple;
    for (auto &Emult  : eigenvalues) {
        tie(E, multiplicity) = Emult;
        eigenvalues_simple.push_back(E);
        const pair<double, double> &error = p2.matchingError(E);
        CHECK(abs(error.first) < 1e-3);

        if (E < 6) {
            double El = p2.eigenvalue(E - 0.01).first;
            double Em = p2.eigenvalue(E + 0.01).first;
            CHECK(Approx(El).margin(1e-7) == E);
            CHECK(Approx(Em).margin(1e-7) == E);
            CHECK(p2.eigenvalueError(El) < 1e-5);
            CHECK(p2.eigenvalueError(Em) < 1e-5);
        }

        const vector<Eigenfunction2D<>> f = p2.eigenfunction(E);
        CHECK(f.size() == multiplicity);
    }
    checkOrthonormality(&p2, eigenvalues_simple.begin(), eigenvalues_simple.end());
}
