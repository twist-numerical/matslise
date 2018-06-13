#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "matscs.h"
#include "matslise.h"
#include "se2d.h"

using namespace std;
using namespace Eigen;

void coffey() {
    double M = M_PI_2;
    double B = 20;

    Matslise coffey([B](double x) {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, -M, M, 128);

    cout << endl;

    matslise::Y y0({0, 1});
    ArrayXd x(6);
    x << 0.0700000000000000,
            0.0760000000000000,
            0.0820000000000000,
            0.0880000000000000,
            0.0940000000000000,
            0.100000000000000;

    vector<tuple<unsigned int, double>> *eigenvalues = coffey.computeEigenvaluesByIndex(0, 2, y0, y0);
    cout.precision(17);
    unsigned int i;
    double E;
    for (tuple<unsigned int, double> iE : *eigenvalues) {
        tie(i, E)  = iE;

        Array<matslise::Y, Dynamic, 1> y = coffey.computeEigenfunction(E, y0, y0, x);
        for (int j=0; j < x.size(); ++j) {
            cout << x[j] << ", ";
        }
        cout << endl;
        for (int j=0; j < y.size(); ++j) {
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

    matslise::Y y0({0, 1});
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

    vector<tuple<unsigned int, double>> *eigenvalues = mathieu.computeEigenvalues(0, 100, y0, y0);
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
    SE2D se2d([](double x, double y) { return (1+x*x)*(1+y*y); }, -5.5,5.5, -5.5,5.5, 30, 15);

    cout << se2d.xmin << endl;
}

int main() {
    coffey();
//    test2d();
}