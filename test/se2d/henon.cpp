#include <cmath>
#include <vector>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "./checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;

const vector<tuple<Index, double, Index>> HENON_EIGENVALUES{
        {0,  2 * 0.998594690530479, 1},
        {1,  2 * 1.99007660445524,  2},
        {3,  2 * 2.95624333869018,  1},
        {4,  2 * 2.98532593386986,  2},
        {6,  2 * 3.92596412795287,  2},
        {8,  2 * 3.98241882458866,  1},
        {9,  2 * 3.98575763690663,  1},
        {10, 2 * 4.87014557482289,  1},
        {11, 2 * 4.89864497284387,  2}};

TEST_CASE("Eigenfunctions henon", "[matslise2d][eigenfunctions][henon][auto]") {
    Matslise2D<>::Config config;
    config.basisSize = 14;
    config.tolerance = 1e-9;
    config.stepsPerSector = 1;

    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * x * (y * y - x * x / 3);
            },
            {-6., 6., -6., 6.}, config);

    checkProblem(problem, HENON_EIGENVALUES);
}

TEST_CASE("Eigenfunctions henon (flipped)", "[matslise2d][eigenfunctions][henon][auto]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-8;

    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * y * (x * x - y * y / 3);
            },
            {-6., 6., -6., 6.}, config);

    checkProblem(problem, HENON_EIGENVALUES);
}

TEST_CASE("Eigenfunctions henon (symmetric)", "[matslise2d][eigenfunctions][henon][symmetric][auto]") {
    Matslise2D<>::Config config;
    config.tolerance = 1e-10;
    config.xSymmetric = true;

    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * y * (x * x - y * y / 3);
            },
            {-6., 6., -6., 6.}, config);

    checkProblem(problem, HENON_EIGENVALUES);
}

TEST_CASE("Eigenfunctions henon (half)", "[matslise2d][eigenfunctions][henon][half][auto]") {
    Matslise2D<>::Config config;
    config.basisSize = 14;
    config.tolerance = 1e-10;
    config.stepsPerSector = 1;

    Matslise2DHalf<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * x * (y * y - x * x / 3);
            },
            {-6., 6., -6., 6.}, config);

    checkProblem(problem, HENON_EIGENVALUES);
}

//https://www.sciencedirect.com/science/article/pii/S0010465508003998
const vector<tuple<Index, double, Index>> HENON_MORE_EIGENVALUES{
        {0,  2 * 0.99859477260462, 1},
        {1,  2 * 1.99007676008316, 2},
        {3,  2 * 2.95624298898781, 1},
        {4,  2 * 2.98532642806441, 2},
        {6,  2 * 3.92596372109113, 2},
        {8,  2 * 3.98241728327456, 1},
        {9,  2 * 3.98576092607719, 1},
        {10, 2 * 4.87014400547251, 1},
        {11, 2 * 4.89864420444642, 2},
        {13, 2 * 4.98625101494171, 2},
        {15, 2 * 5.81701909971053, 2},
        {17, 2 * 5.86701480916528, 1},
        {18, 2 * 5.88144609873667, 1},
        {19, 2 * 5.99132695571450, 2},
        {21, 2 * 6.73791623075037, 1},
        {22, 2 * 6.76486656303911, 2},
        {24, 2 * 6.85343062732831, 2},
        {26, 2 * 6.99893192820090, 1},
        {27, 2 * 6.99938690825129, 1},
        {28, 2 * 7.65948550690946, 2}
        /*  {30, 2 * 7.69772136553152, 1},
          {31, 2 * 7.73688473693607, 1},
          {32, 2 * 7.83273518682762, 2},
          {34, 2 * 8.00942477463143, 2},
          {36, 2 * 8.55402322299586, 1},
          {37, 2 * 8.57635148658287, 2},
          {39, 2 * 8.67792887113521, 2},
          {41, 2 * 8.81132713081544, 1},
          {42, 2 * 8.81518847088363, 1},
          {43, 2 * 9.02172330707064, 2},
          {45, 2 * 9.44405461567860, 1},*/
};


TEST_CASE("Eigenfunctions henon extended", "[matslise2d][eigenfunctions][henon][half][auto][slow]") {
    Matslise2D<>::Config config;
    config.basisSize = 30;
    config.tolerance = 1e-9;
    config.stepsPerSector = 1;
    config.xSymmetric = true;

    Matslise2D<> problem(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2. * sqrt(5.)) * y * (x * x - y * y / 3);
            },
            {-10., 10., -10., 10.}, config);

    checkProblem(problem, HENON_MORE_EIGENVALUES);
}