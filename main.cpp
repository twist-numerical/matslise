#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <matscs.h>
#include "matslise/matslise.h"

using namespace std;
using namespace Eigen;

int main() {
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

    vector<double> *eigenvalues = mathieu.computeEigenvalues(0, 100, y0, y0);
    cout.precision(17);
    for (double E : *eigenvalues) {
        cout << E << endl;
    }
    delete eigenvalues;

    Matslise coffey([](double x) {
        return -2 * 30 * cos(2 * x) + 30 * 30 * sin(2 * x) * sin(2 * x);
    }, m, M, 256);

    eigenvalues = coffey.computeEigenvalues(0, 100, y0, y0);
    cout.precision(17);
    for (double E : *eigenvalues) {
        cout << E << endl;
    }
    delete eigenvalues;
    return 0;
}