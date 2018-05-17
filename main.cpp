#include <iostream>
#include <cmath>
#include <vector>
#include <matslise.h>

using namespace std;

int main() {
    Matslise ms(matslise::Mathieu(), 0, M_PI, 1000);
    double E = 400.001253135;
    vector<double> x;
    int n = 10001;
    for (int i = 0; i < n; ++i)
        x.push_back(i * M_PI / (n - 1));
    vector<matslise::Y> *ys = ms.computeEigenfunction(E, x);
    int i = 0;
    for (auto &v  : *ys) {
        if (++i % 100 == 0)
            cout << v << ", ";
    }

    delete ys;

    return 0;
}