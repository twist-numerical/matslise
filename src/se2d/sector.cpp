//
// Created by toon on 6/13/18.
//

#include <iostream>
#include "../se2d.h"
#include "../util/lobatto.h"

using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace matslise::se2d_util;

template<typename FROM, typename TO>
Array<TO, Dynamic, 1> apply(const Array<FROM, Dynamic, 1> &x, function<TO(FROM)> f) {
    Array<TO, Dynamic, 1> r(x.size());
    for (int i = 0; i < r.size(); ++i)
        r[i] = f(x[i]);
    return r;
}


Sector::Sector(SE2D *se2d, double ymin, double ymax, int matslise_count, int matscs_count)
        : se2d(se2d), ymin(ymin), ymax(ymax) {
    const Y<double> y0 = Y<double>({0, 1});

    const double ybar = (ymax + ymin) / 2;
    function<double(double)> vbar_fun = [se2d, ybar](double x) -> double { return se2d->V(x, ybar); };
    vbar = apply(se2d->xGrid, vbar_fun);
    matslise = new Matslise(vbar_fun, se2d->xmin, se2d->xmax, matslise_count);

    vector<pair<unsigned int, double>> *index_eigv = matslise->computeEigenvaluesByIndex(0, se2d->N, y0, y0);
    eigenvalues = new double[se2d->N];
    eigenfunctions = new ArrayXd[se2d->N];
    eigenfunctionsScaling = new double[se2d->N];
    for (int i = 0; i < se2d->N; ++i) {
        double E = (*index_eigv)[i].second;
        eigenvalues[i] = E;
        Array<Y<double>, Dynamic, 1> func = matslise->computeEigenfunction(E, y0, y0, se2d->xGrid);
        eigenfunctions[i] = ArrayXd(func.size());
        for (int j = 0; j < func.size(); ++j)
            eigenfunctions[i][j] = func[j].y[0];
        eigenfunctions[i] *= (eigenfunctionsScaling[i] = 1 /
                sqrt(lobatto::quadrature(se2d->xGrid, eigenfunctions[i] * eigenfunctions[i])));
    }

    matscs = new Matscs([this](double y) -> MatrixXd { return this->calculateDeltaV(y); }, se2d->N, ymin, ymax,
                        matscs_count);

    delete index_eigv;
}

Sector::~Sector() {
    delete matslise;
    delete matscs;
    delete[] eigenvalues;
    delete[] eigenfunctions;
    delete[] eigenfunctionsScaling;
}

Y<MatrixXd> Sector::propagate(double E, const Y<MatrixXd> &c, bool forward) const {
    return propagate(E, c, forward ? ymax : ymin, forward);
}

Y<MatrixXd> Sector::propagate(double E, const Y<MatrixXd> &c, double y, bool forward) const {
    return matscs->propagate(E, c, forward ? ymin : ymax, y);
}

MatrixXd Sector::calculateDeltaV(double y) const {
    Eigen::MatrixXd dV(se2d->N, se2d->N);

    ArrayXd vDiff =
            apply<double, double>(se2d->xGrid, [this, y](double x) -> double { return this->se2d->V(x, y); }) - vbar;

    for (int i = 0; i < se2d->N; ++i) {
        for (int j = 0; j <= i; ++j) {
            dV(i, j) = lobatto::quadrature(se2d->xGrid, eigenfunctions[i] * eigenfunctions[j] * vDiff);
            if (j < i) dV(j, i) = dV(i, j);
        }
        dV(i, i) += eigenvalues[i];
    }

    return dV;
}

ArrayXd Sector::computeEigenfunction(int index, const ArrayXd &x) const {
    const Y<double> y0 = Y<double>({0, 1});
    long n = x.size();

    Array<matslise::Y<double>, Dynamic, 1> raw = matslise->computeEigenfunction(eigenvalues[index], y0, y0, x);
    ArrayXd result(n);
    for (int i = 0; i < n; ++i)
        result(i) = raw(i).y[0] * eigenfunctionsScaling[index];
    return result;
}