//
// Created by toon on 6/13/18.
//

#include <iostream>
#include <map>
#include "../se2d.h"
#include "../util/lobatto.h"

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
SEBase<n>::SEBase(typename dim<n>::function V, const matslise::Rectangle<n> &domain, const Options<n> &options) :
        V(V), domain(domain),
        sectorCount(options._sectorCount), N(options._N) {
    sectors = new Sector<n> *[sectorCount];

    for (int i = 0; i < n - 1; ++i) {
        grid[i] = lobatto::grid(getGrid(domain.getMin(i), domain.getMax(i), options._gridPoints));
    }

    double h = (domain.max - domain.min) / sectorCount;
    for (int i = 0; i < sectorCount; ++i)
        sectors[i] = new Sector<n>(this, domain.min + i * h, domain.min + (i + 1) * h, options);

    M = new MatrixXd[sectorCount - 1];
    for (int i = 0; i < sectorCount - 1; ++i)
        M[i] = calculateM(i);
}

template<int n>
MatrixXd SEBase<n>::calculateM(int k) const {
    MatrixXd M(N, N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            M(i, j) = lobatto::multi_quadrature<n - 1>(grid, sectors[k]->eigenfunctions[j] *
                                                             sectors[k + 1]->eigenfunctions[i]);

    // MatrixXd Q = M.householderQr().householderQ();
    return M;
}


template<int n>
SEBase<n>::~SEBase() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
    delete[] M;
}

template<int n>
pair<vector<MatrixXd>, vector<MatrixXd>> SEBase<n>::calculateAllSteps(double E) const {
    vector<MatrixXd> y(sectorCount - 1);
    vector<MatrixXd> dy(sectorCount - 1);

    int match = sectorCount / 2;
    Y<MatrixXd> y0({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)}, {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});
    Y<MatrixXd> yl = y0;
    for (int i = 0; i <= match; ++i) {
        yl = M[i] * sectors[i]->propagate(E, yl, true);
        y[i] = yl.y[0];
        dy[i] = yl.y[1];
        yl *= (yl.y[0].colwise().norm()).cwiseInverse().asDiagonal();
    }
    Y<MatrixXd> yr = sectors[sectorCount - 1]->propagate(E, y0, false);
    for (int i = sectorCount - 2; i > match; --i) {
        yr = sectors[i]->propagate(E, M[i].transpose() * yr, false);
        y[i] = yr.y[0];
        dy[i] = yr.y[1];
        yr *= (yr.y[0].colwise().norm()).cwiseInverse().asDiagonal();
    }
    return make_pair(y, dy);
}

template<int n>
pair<MatrixXd, MatrixXd> SEBase<n>::calculateErrorMatrix(double E) const {
    int match = sectorCount / 2;
    Y<MatrixXd> y0({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)}, {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});
    Y<MatrixXd> yl = y0;
    for (int i = 0; i <= match; ++i) {
        yl = M[i] * sectors[i]->propagate(E, yl, true);
        yl *= (yl.y[0].colwise().norm()).cwiseInverse().asDiagonal();
    }
    Y<MatrixXd> yr = sectors[sectorCount - 1]->propagate(E, y0, false);
    for (int i = sectorCount - 2; i > match; --i) {
        yr = sectors[i]->propagate(E, M[i].transpose() * yr, false);
        yr *= (yr.y[0].colwise().norm()).cwiseInverse().asDiagonal();
    }


    ColPivHouseholderQR<MatrixXd> left_solver(yl.y[0].transpose());
    ColPivHouseholderQR<MatrixXd> right_solver(yr.y[0].transpose());
    MatrixXd Ul = left_solver.solve(yl.y[1].transpose()).transpose();
    MatrixXd Ur = right_solver.solve(yr.y[1].transpose()).transpose();
    return make_pair(
            Ul - Ur,
            left_solver.solve((yl.dy[1] - Ul * yl.dy[0]).transpose()).transpose()
            - right_solver.solve((yr.dy[1] - Ur * yr.dy[0]).transpose()).transpose()
    );
}

template<int n>
vector<pair<double, double>> *SEBase<n>::calculateErrors(double E) const {
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
SEBase<n>::sortedErrors(double E,
                        const std::function<bool(std::pair<double, double>, std::pair<double, double>)> &sorter) const {
    vector<pair<double, double>> *errors = calculateErrors(E);

    sort(errors->begin(), errors->end(), sorter);

    return errors;
}

template<int n>
double SEBase<n>::findEigenvalue(double E) const {
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
pair<double, double> SEBase<n>::calculateError(
        double E, const std::function<bool(std::pair<double, double>, std::pair<double, double>)> &sorter) const {
    vector<pair<double, double>> *errors = sortedErrors(E, sorter);
    pair<double, double> best = (*errors)[0];
    delete errors;
    return best;
}

template
class matslise::SEBase<2>;

template
class matslise::SEBase<3>;