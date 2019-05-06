//
// Created by toon on 6/13/18.
//

#include <iostream>
#include <map>
#include "../se2d.h"
#include "../util/lobatto.h"
#include "../util/find_sector.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;


ArrayXd getGrid(double min, double max, int count) {
    ArrayXd points(count);
    for (int i = 0; i < count; ++i)
        points[i] = min + (max - min) * i / (count - 1);
    return points;
}

template<int n>
SEnD<n>::SEnD(typename dim<n>::function V, const matslise::Rectangle<n> &domain, const Options<n> &_options) :
        V(V), domain(domain), N(_options._N),
        options(_options) {
    for (int i = 0; i < n - 1; ++i) {
        grid[i] = lobatto::grid(getGrid(domain.getMin(i), domain.getMax(i), options._gridPoints));
    }
    options._builder->build(this, domain.min, domain.max);

    M = new MatrixXd[sectorCount - 1];
    for (int i = 0; i < sectorCount - 1; ++i)
        M[i] = calculateM(i);
}

template<int n>
MatrixXd SEnD<n>::calculateM(int k) const {
    MatrixXd M(N, N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            M(i, j) = lobatto::multi_quadrature<n - 1>(grid, sectors[k]->eigenfunctions[j] *
                                                             sectors[k + 1]->eigenfunctions[i]);

    // MatrixXd Q = M.householderQr().householderQ();
    return M;
}


template<int n>
SEnD<n>::~SEnD() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
    delete[] M;
}

template<int n>
pair<MatrixXd, MatrixXd> SEnD<n>::calculateErrorMatrix(double E) const {
    Y<Dynamic> y0 = Y<Dynamic>::Dirichlet(N);
    Y<Dynamic> yl = y0;
    for (int i = 0; sectors[i]->max <= match; ++i) {
        yl = M[i] * sectors[i]->propagate(E, yl, sectors[i]->min, sectors[i]->max, true);
        //yl *= (yl.getY(0).colwise().norm()).cwiseInverse().asDiagonal();
    }
    Y<Dynamic> yr = sectors[sectorCount - 1]->propagate(
            E, y0, sectors[sectorCount - 1]->max, sectors[sectorCount - 1]->min, true);
    for (int i = sectorCount - 2; sectors[i]->min >= match; --i) {
        yr = sectors[i]->propagate(E, (MatrixXd) (M[i].transpose()) * yr, sectors[i]->max, sectors[i]->min, true);
        //yr *= (yr.getY(0).colwise().norm()).cwiseInverse().asDiagonal();
    }


    ColPivHouseholderQR<MatrixXd> left_solver(yl.getY(0).transpose());
    ColPivHouseholderQR<MatrixXd> right_solver(yr.getY(0).transpose());
    MatrixXd Ul = left_solver.solve(yl.getY(1).transpose()).transpose();
    MatrixXd Ur = right_solver.solve(yr.getY(1).transpose()).transpose();
    return make_pair(
            Ul - Ur,
            left_solver.solve((yl.getdY(1) - Ul * yl.getdY(0)).transpose()).transpose()
            - right_solver.solve((yr.getdY(1) - Ur * yr.getdY(0)).transpose()).transpose()
    );
}

template<int n>
Y<Dynamic> SEnD<n>::propagate(double E, const Y<Dynamic> &y0, double a, double b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Dynamic> y = y0;
    int sectorIndex = find_sector<SEnD<n>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        if (direction == -1 && sectorIndex < sectorCount - 1)
            y = (MatrixXd) (M[sectorIndex].transpose()) * y;
        y = sector->propagate(E, y, a, b, use_h);
        if (direction == 1)
            y = M[sectorIndex] * y;
        sectorIndex += direction;
    } while (!sector->contains(b));
    return y;
}

template<int n>
vector<pair<double, double>> *SEnD<n>::calculateErrors(double E) const {
    pair<MatrixXd, MatrixXd> error_matrix = calculateErrorMatrix(E);
    EigenSolver<MatrixXd> solver(N);

    solver.compute(error_matrix.first, true);

    multimap<double, int> rightMap;
    ArrayXcd eigenvaluesRight = solver.eigenvalues().array();
    for (int i = 0; i < N; ++i)
        rightMap.insert({eigenvaluesRight[i].real(), i});
    MatrixXcd right = solver.eigenvectors();

    multimap<double, int> leftMap;
    solver.compute(error_matrix.first.transpose(), true);
    ArrayXcd eigenvaluesLeft = solver.eigenvalues().array();
    for (int i = 0; i < N; ++i)
        leftMap.insert({eigenvaluesLeft[i].real(), i});
    MatrixXcd left = solver.eigenvectors().transpose();


    auto errors = new vector<pair<double, double>>();
    for (auto leftI = leftMap.begin(), rightI = rightMap.begin();
         leftI != leftMap.end() && rightI != rightMap.end();
         ++leftI, ++rightI) {
        int &li = leftI->second;
        int &ri = rightI->second;
        errors->push_back({eigenvaluesLeft[li].real(),
                           (left.row(li) * error_matrix.second * right.col(ri) /
                            (left.row(li) * right.col(ri)))[0].real()});
    }

    return errors;
}

template<int n>
vector<pair<double, double>> *
SEnD<n>::sortedErrors(double E,
                      const std::function<bool(std::pair<double, double>, std::pair<double, double>)> &sorter) const {
    vector<pair<double, double>> *errors = calculateErrors(E);

    sort(errors->begin(), errors->end(), sorter);

    return errors;
}

template<int n>
double SEnD<n>::findEigenvalue(double E) const {
    double error, derror;
    int i = 0;
    do {
        tie(error, derror) = calculateError(E, NEWTON_RAPHSON_SORTER);
        E -= error / derror;
        ++i;
    } while (i < 30 && abs(error) > 1e-9);

    if (abs(error) > 1e-5)
        return NAN;
    return E;
}

template<int n>
pair<double, double> SEnD<n>::calculateError(
        double E, const std::function<bool(std::pair<double, double>, std::pair<double, double>)> &sorter) const {
    vector<pair<double, double>> *errors = sortedErrors(E, sorter);
    pair<double, double> best = (*errors)[0];
    delete errors;
    return best;
}

template
class matslise::SEnD<2>;

template
class matslise::SEnD<3>;