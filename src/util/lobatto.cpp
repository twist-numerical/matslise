#include "lobatto.h"

using namespace Eigen;

ArrayXd lobatto::grid(const ArrayXd &x) {
    long n = x.size();
    ArrayXd fine(3 * n - 2);

    for (long i = 0; i < n - 1; ++i) {
        double a = x[i], b = x[i + 1];
        double m = (a + b) / 2, h = (b - a) / 2;
        fine[3 * i] = a;
        fine[3 * i + 1] = m - h / sqrt(5);
        fine[3 * i + 2] = m + h / sqrt(5);
    }
    fine[3 * (n - 1)] = x[n - 1];

    return fine;
}

double lobatto::quadrature(const ArrayXd &x, const ArrayXd &f) {
    long n = x.size();
    double result = 0;

    for (int i = 0; i < n - 1; i += 3)
        result += (x[i + 3] - x[i]) / 2 * ((f[i] + f[i + 3]) / 6 + (f[i + 1] + f[i + 2]) * 5 / 6);

    return result;
}

template<>
double lobatto::multi_quadrature<1>(const ArrayXd x[1], const typename Dim<1>::grid &f) {
    return lobatto::quadrature(x[0], f);
}

template<>
double lobatto::multi_quadrature<2>(const ArrayXd x[2], const typename Dim<2>::grid &f) {
    long n = x[0].size();
    ArrayXd result = ArrayXd::Zero(n);

    for (int i = 0; i < n - 1; i += 3)
        result += (x[0][i + 3] - x[0][i]) / 2 * ((f.col(i) + f.col(i + 3)) / 6 + (f.col(i + 1) + f.col(i + 2)) * 5 / 6);

    return lobatto::quadrature(x[1], result);
}