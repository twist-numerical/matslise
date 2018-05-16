#include <iostream>
#include <cmath>
#include "matslise.h"

using namespace std;

int main() {
    Matslise ms(matslise::Mathieu(), -M_PI, M_PI, 16);
    double h = M_PI / 30;
    matslise::Y y = {0, 1};

    cout << ms.sectors[0]->calculateT(3.91702477299847, ms.sectors[0]->h) << endl;
    cout << ms.propagate(3.91702477299847, y, 0, M_PI) << endl;

    cout << "(" << M_PI << ", " << y << ")";
    for (double d = M_PI; d - h / 2 > -M_PI; d -= h) {
        y = ms.propagate(3.91702477299847, y, d, d - h);
        cout << ", (" << (d-h) << ", " << y << ")";
    }
    cout << endl;
    return 0;
}