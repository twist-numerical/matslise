#include <iostream>
#include <cmath>
#include <vector>
#include "matslise/matscs.h"

using namespace std;
using namespace Eigen;

int main() {
    Matscs ms([](double x) -> MatrixXd {
        MatrixXd m(2, 2);
        m << 3*x, -x, -x, 3*x;
        return m;
    }, 2, 0, 1, 16);

    double E = 10.368507161836;

    cout << ms.propagate(E, matscs::Y(MatrixXd::Zero(2,2), MatrixXd::Identity(2,2)), 0, 1) << endl;
    /*
    vector<matscs::Y> *ys = ms.computeEigenfunction(E, x);
    int i = 0;
    for (auto &v  : *ys) {
        if (++i % 100 == 0)
            cout << v << ", ";
    }

    delete ys;
     */

    return 0;
}