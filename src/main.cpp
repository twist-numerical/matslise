#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "matscs.h"
#include "matslise.h"

using namespace std;
using namespace Eigen;

void coffey() {
    double M = M_PI_2;

    matslise::HalfRange coffey([](double x) {
        return -2 * 30 * cos(2 * x) + 30 * 30 * sin(2 * x) * sin(2 * x);
    }, M, 60);

    cout << endl;

    matslise::Y y0({0, 1});
    vector<double> x = {-.6,.1,1.1};
    vector<double> *eigenvalues = coffey.computeEigenvalues(0, 1000, y0);
    cout.precision(17);
    for (double E : *eigenvalues) {
        coffey.computeEigenfunction(E, y0, x);
        cout << E << endl;
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

    vector<double> *eigenvalues = mathieu.computeEigenvalues(0, 100, y0, y0);
    cout.precision(17);
    for (double E : *eigenvalues) {
        cout << E << endl;
    }
    delete eigenvalues;
}


int main() {
    coffey();
}