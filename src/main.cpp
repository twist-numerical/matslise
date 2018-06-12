#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "matscs.h"
#include "matslise.h"

using namespace std;
using namespace Eigen;

void coffey() {
    double M = M_PI_2;
    double B = 20;

    Matslise coffey([B](double x) {
        return B*B*x-2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, -M, M, 256);

    cout << endl;

    matslise::Y y0({0, 1});
    vector<double> x = {-1,1};

    vector<tuple<unsigned int, double>> *eigenvalues = coffey.computeEigenvaluesByIndex(0, 2, y0, y0);
    cout.precision(17);
    unsigned int i;
    double E;
    for (tuple<unsigned int, double> iE : *eigenvalues) {
        tie(i, E)  = iE;

        vector<matslise::Y> *ys = coffey.computeEigenfunction(E, y0, y0, x);
        for (double _x : x) {
            cout << _x << ", ";
        }
        cout << endl;
        for (matslise::Y &y : *ys) {
            cout << y.y[0] << ", ";
        }
        delete ys;

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


int main() {
    coffey();
}