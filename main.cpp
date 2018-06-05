#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <matscs.h>
#include "matslise/matslise.h"

using namespace std;
using namespace Eigen;

int main() {
    double m = -M_PI_2, M = M_PI_2;

    double E = 10;
    Matslise coffey([](double x) {
        return 2*cos(2*x);
        //return -2 * 30 * cos(2 * x) + 30 * 30 * sin(2 * x) * sin(2 * x);
    }, m, M, 16);

    matslise::Y y1, y0 = matslise::Y(0, 1);
    double theta;
    tie(y1, theta) = coffey.propagate(E, y0, m, M);
    cout << theta << endl;
    cout << coffey.calculateError(E, y0, y0) << endl;

    return 1;
    Matscs ms([](double x) -> MatrixXd {
        MatrixXd m(2, 2);
        m << 3 * x, -x, -x, 3 * x;
        return m;
    }, 2, 0, 1, 16);

    double matE = 10.368507161836;
    const matscs::Y &y = matscs::Y(MatrixXd::Zero(2, 2), MatrixXd::Identity(2, 2));
    cout << ms.propagate(matE, y, 0, 1) << endl;


    vector<double> x;
    int n = 500;
    for (int i = 0; i <= 100; ++i)
        x.push_back(M_PI + 2 * M_PI * i / n);

    vector<matslise::Y> *ys = coffey.computeEigenfunction(E, x);
    int i = 0;
    for (auto &v  : *ys) {
        if (++i % 10 == 0)
            cout << v << ", ";
    }

    delete ys;


    return 0;
}