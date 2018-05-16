#include <iostream>
#include <cmath>
#include "matslise/matslise.h"

double Mathieu(double x) {
    return 2 * cos(2 * x);
}

using namespace std;

int main() {
    Matslise ms(&Mathieu, 0, M_PI, 16);
    double h = M_PI / 30;
    Matslise::Y y = {0, 1};

    cout << ms.sectors[0]->calculateT(3.91702477299847, ms.sectors[0]->h) << endl;
    cout << ms.propagate(3.91702477299847, y, 0, M_PI) << endl;

    std::cout << "(" << 0 << ", " << y << ")";
    for (double d = 0; d + h / 2 < M_PI; d += h) {
        y = ms.propagate(3.91702477299847, y, d, d + h);
        std::cout << ", (" << (d + h) << ", " << y << ")";
    }
    std::cout << std::endl;
    return 0;
}