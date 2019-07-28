#include <iostream>
#include <cmath>
#include "../src/matslise.h"

using namespace std;
using namespace Eigen;
using namespace matslise;

void test2d() {
    std::cout.precision(10);
    SE2D<> se2d([](double x, double y) { return (1 + x * x) * (1 + y * y); }, {{-5.5, 5.5}, -5.5, 5.5},
                Options2<>().sectorCount(11).nested(Options1<>().sectorCount(26)).N(6));
    cout << se2d.propagate(3, Y<double, -1>::Dirichlet(se2d.N), -5.5, 0).getY(0) << endl;
}

int main() {
    test2d();
}