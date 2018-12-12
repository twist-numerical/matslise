#include <iostream>
#include <cmath>
#include <matslise/matscs.h>
#include <matslise/se2d.h>

using namespace std;
using namespace Eigen;
using namespace matslise;

void coffey() {
    double M = M_PI_2;
    double B = 20;

    matslise::Matslise coffey([B](double x) {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, -M, M, 128);

    cout << endl;

    matslise::Y<double> y0({0, 1});
    ArrayXd x(6);
    x << 0.0700000000000000,
            0.0760000000000000,
            0.0820000000000000,
            0.0880000000000000,
            0.0940000000000000,
            0.100000000000000;

    vector<pair<int, double>> *eigenvalues = coffey.computeEigenvaluesByIndex(0, 2, y0, y0);
    cout.precision(17);
    unsigned int i;
    double E;
    for (pair<int, double> iE : *eigenvalues) {
        tie(i, E) = iE;

        Array<matslise::Y<double>, Dynamic, 1> y = coffey.computeEigenfunction(E, y0, y0, x);
        for (int j = 0; j < x.size(); ++j) {
            cout << x[j] << ", ";
        }
        cout << endl;
        for (int j = 0; j < y.size(); ++j) {
            cout << y[j].y[0] << ", ";
        }

        cout << endl << i << ": " << E << endl << endl;
    }
    delete eigenvalues;
}

void mathieu() {
    double m = -M_PI_2, M = M_PI_2;

    Matslise mathieu([](double x) {
        return 2 * cos(2 * x);
        //return -2 * 30 * cos(2 * x) + 30 * 30 * sin(2 * x) * sin(2 * x);
    }, m, M, 16);
    double h = M_PI / 16;

    matslise::Y<double> y0({0, 1});
    cout << get<0>(mathieu.propagate(110, y0, m, m + 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(120, y0, m, m + 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(130, y0, m, m + 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(140, y0, m, m + 2 * h)) << endl;
    cout << endl;
    cout << get<0>(mathieu.propagate(110, y0, M, M - 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(120, y0, M, M - 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(130, y0, M, M - 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(140, y0, M, M - 2 * h)) << endl;
    cout << endl;

    vector<pair<int, double>> *eigenvalues = mathieu.computeEigenvalues(0, 100, y0, y0);
    cout.precision(17);

    unsigned int i;
    double E;
    for (auto &iE : *eigenvalues) {
        tie(i, E) = iE;
        cout << i << ": " << E << endl;
    }
    delete eigenvalues;
}


void test2d() {
    std::cout.precision(10);
    SEnD<2> se2d([](double x, double y) { return (1 + x * x) * (1 + y * y); }, {{-5.5, 5.5}, -5.5, 5.5},
                 Options<2>().sectorCount(64).nested(Options<1>().sectorCount(64)).N(8));

    /*for(double a = 3.1; a < 3.3; a += .001)
        cout << se2d.calculateError(a).first << endl;*/
    ArrayXd xs(5);
    xs << -2, -1, 0, 1, 2;
    ArrayXd ys(5);
    ys << -2, -1, 0, 1, 2;
    for (auto &a : se2d.computeEigenfunction(5.535959903292571, xs, ys))
        cout << a << "\n\n" << endl;
}


void test3d() {
    std::cout.precision(10);
    SEnD<3> se3d([](double x, double y, double z) -> double { return (1 + x * x) * (1 + y * y) * (1 + z * z); },
                 {{{-5.5, 5.5}, -5.5, 5.5}, -5.5, 5.5},
                 Options<3>().sectorCount(25).stepsPerSector(5).N(8).nested(
                         Options<2>().sectorCount(25).stepsPerSector(5).N(8).nested(
                                 Options<1>().sectorCount(35))));

    for (double a = 3.1; a < 3.3; a += .01)
        cout << se3d.calculateError(a).first << endl;
}

void testMatscs() {
    Matscs ms([](double x) {
        return (MatrixXd(1, 1) << x * x).finished();
    }, 1, -5, 5, 32);
    MatrixXd zero = MatrixXd::Zero(1, 1);
    MatrixXd one = MatrixXd::Identity(1, 1);
    double E = 3;
    cout << ms.propagate(E, Y<MatrixXd>({zero, one}, {zero, zero}), -5, 5).y[0] << endl;

    cout << "\n\nMATSLISE" << endl;
    Matslise ml([](double x) { return x * x; }, -5, 5, 32);
    cout << get<0>(ml.propagate(E, Y<double>({0, 1}, {0, 0}), -5, 5)).y[0] << endl;
}

void testEigenfunctionCalculator() {
    double M = M_PI_2;
    double B = 20;

    matslise::Matslise coffey([B](double x) {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, -M, M, 128);

    Y<double> y0({0, 1});

    vector<pair<int, double>> *eigs = coffey.computeEigenvaluesByIndex(2, 3, y0, y0);
    double E = get<1>((*eigs)[0]);
    cout << E << endl;
    delete eigs;

    matslise_util::EigenfunctionCalculator *_ef = coffey.eigenfunctionCalculator(E, y0, y0);
    matslise_util::EigenfunctionCalculator &ef = *_ef;
    cout << ef(0) << endl;
    delete _ef;
}

void testPrufer() {
    //Matslise ms([](double x) { return 0.0; }, 0., 3.14, 60);

    Matslise ms([](double x) {
        return 2 * 5 * cos(2 * x);
    }, -M_PI_2, M_PI_2, 18);
    matslise::Y<double> y0({0, 1});
    cout << get<1>(ms.propagate(1, y0, -M_PI_2, -1)) << endl;
    // cout << get<0>(ms.propagate(16, y0, .1, 1.6)) << endl;
}

int main() {
    // coffey();
    // testPrufer();
    test3d();
    // testBigE();
    // testMatscs();
    // testEigenfunctionCalculator();
}