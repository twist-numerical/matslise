//
// Created by toon on 6/13/18.
//

#include <iostream>
#include "../se2d.h"
#include "../util/lobatto.h"

using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace matslise::SEnD_util;

template<int n>
typename dim<n>::array apply(const ArrayXd grid[n], const typename dim<n>::function &f);

template<>
typename dim<1>::array apply<1>(const ArrayXd grid[1], const dim<1>::function &f) {
    ArrayXd r(grid[0].size());
    for (int i = 0; i < r.size(); ++i)
        r[i] = f(grid[0][i]);
    return r;
}

template<>
typename dim<2>::array apply<2>(const ArrayXd grid[2], const dim<2>::function &f) {
    ArrayXXd r(grid[1].size(), grid[0].size());
    for (int i = 0; i < grid[1].size(); ++i)
        for (int j = 0; j < grid[0].size(); ++j)
            r(i, j) = f(grid[0][j], grid[1][i]);
    return r;
}

template<>
Sector<2>::Sector(SEBase<2> *se2d, double ymin, double ymax, const Options<2> &options)
        : se2d(se2d), min(ymin), max(ymax) {
    const Y<double> y0 = Y<double>({0, 1}, { 0,0 });

    const double ybar = (ymax + ymin) / 2;
    function<double(double)> vbar_fun = [se2d, ybar](double x) -> double { return se2d->V(x, ybar); };
    vbar = apply<1>(se2d->grid, vbar_fun);
    matslise = new Matslise(vbar_fun, se2d->domain.sub.min, se2d->domain.sub.max, options.nestedOptions._sectorCount);

    vector<pair<int, double>> *index_eigv = matslise->computeEigenvaluesByIndex(0, se2d->N, y0, y0);
    eigenvalues = new double[se2d->N];
    eigenfunctions = new ArrayXd[se2d->N];
    eigenfunctionsScaling = new double[se2d->N];
    for (int i = 0; i < se2d->N; ++i) {
        double E = (*index_eigv)[i].second;
        eigenvalues[i] = E;
        Array<Y<double>, Dynamic, 1> func = matslise->computeEigenfunction(E, y0, y0, se2d->grid[0]);
        eigenfunctions[i] = ArrayXd(func.size());
        for (int j = 0; j < func.size(); ++j)
            eigenfunctions[i][j] = func[j].y[0];
        eigenfunctions[i] *= (eigenfunctionsScaling[i] = 1 /
                                                         sqrt(lobatto::quadrature(se2d->grid[0], eigenfunctions[i] *
                                                                                                 eigenfunctions[i])));
    }

    matscs = new Matscs([this](double y) -> MatrixXd { return this->calculateDeltaV(y); }, se2d->N, ymin, ymax,
                        options._stepsPerSector);

    delete index_eigv;
}

template<>
Sector<3>::Sector(SEBase<3> *se2d, double zmin, double zmax, const Options<3> &options)
        : se2d(se2d), min(zmin), max(zmax) {
    const double zbar = (zmax + zmin) / 2;
    function<double(double, double)> vbar_fun = [se2d, zbar](double x, double y) -> double {
        return se2d->V(x, y, zbar);
    };
    vbar = apply<2>(se2d->grid, vbar_fun);
    matslise = new SEnD<2>(vbar_fun, se2d->domain.sub, options.nestedOptions);

    vector<double> *index_eigv = matslise->computeEigenvaluesByIndex(0, se2d->N);
    eigenvalues = new double[se2d->N];
    eigenfunctions = new ArrayXXd[se2d->N];
    eigenfunctionsScaling = new double[se2d->N];
    for (int i = 0; i < se2d->N;) {
        double E = (*index_eigv)[i];
        eigenvalues[i] = E;
        std::vector<typename dim<2>::array>* funcs = matslise->computeEigenfunction(E, se2d->grid[0], se2d->grid[1]);
        for (auto func : *funcs) {
            eigenfunctions[i] = func;
            eigenfunctions[i] *= (
                    eigenfunctionsScaling[i] =
                            1. / sqrt(lobatto::multi_quadrature<2>(se2d->grid, eigenfunctions[i] * eigenfunctions[i]))
            );
            if (++i >= se2d->N)
                break;
        }
        delete funcs;
    }

    delete index_eigv;

    matscs = new Matscs([this](double y) -> MatrixXd { return this->calculateDeltaV(y); }, se2d->N, zmin, zmax,
                        options._stepsPerSector);

}

template<int n>
Sector<n>::~Sector() {
    delete matslise;
    delete matscs;
    delete[] eigenvalues;
    delete[] eigenfunctions;
    delete[] eigenfunctionsScaling;
}

template<int n>
Y<MatrixXd> Sector<n>::propagate(double E, const Y<MatrixXd> &c, bool forward) const {
    return propagate(E, c, forward ? max : min, forward);
}

template<int n>
Y<MatrixXd> Sector<n>::propagate(double E, const Y<MatrixXd> &c, double y, bool forward) const {
    return matscs->propagate(E, c, forward ? min : max, y);
}

template<int n>
typename dim<n - 1>::function fillIn(typename dim<n>::function &, double);

template<>
typename dim<1>::function fillIn<2>(typename dim<2>::function &f, double y) {
    return [f, y](double x) -> double { return f(x, y); };
}

template<>
typename dim<2>::function fillIn<3>(typename dim<3>::function &f, double z) {
    return [f, z](double x, double y) -> double { return f(x, y, z); };
}

template<int n>
MatrixXd Sector<n>::calculateDeltaV(double z) const {
    Eigen::MatrixXd dV(se2d->N, se2d->N);

    typename dim<n - 1>::array vDiff = apply<n - 1>(se2d->grid, fillIn<n>(this->se2d->V, z)) - vbar;

    for (int i = 0; i < se2d->N; ++i) {
        for (int j = 0; j <= i; ++j) {
            dV(i, j) = lobatto::multi_quadrature<n - 1>(se2d->grid, eigenfunctions[i] * eigenfunctions[j] * vDiff);
            if (j < i) dV(j, i) = dV(i, j);
        }
        dV(i, i) += eigenvalues[i];
    }

    return dV;
}

template<>
ArrayXd Sector<2>::computeEigenfunction(int index, const ArrayXd &x) const {
    const Y<double> y0 = Y<double>({0, 1}, { 0,0 });
    long size = x.size();

    Array<matslise::Y<double>, Dynamic, 1> raw = matslise->computeEigenfunction(eigenvalues[index], y0, y0, x);
    ArrayXd result(size);
    for (int i = 0; i < size; ++i)
        result(i) = raw(i).y[0] * eigenfunctionsScaling[index];
    return result;
}

template
class matslise::SEnD_util::Sector<2>;

template
class matslise::SEnD_util::Sector<3>;