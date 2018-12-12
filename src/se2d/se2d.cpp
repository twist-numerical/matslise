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
SEBase<n>::SEBase(typename dim<n>::function V, const Rectangle<n> &domain, const Options<n> &options) :
        V(V), domain(domain),
        sectorCount(options._sectorCount), N(options._N) {
    sectors = new Sector<n> *[sectorCount];

    for(int i = 0; i < n-1; ++i) {
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
            M(i, j) = lobatto::quadrature(grid[0], sectors[k]->eigenfunctions[j] * sectors[k + 1]->eigenfunctions[i]);

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
pair<MatrixXd, MatrixXd> SEBase<n>::calculateErrorMatrix(double E) const {
    int match = sectorCount / 2;
    Y<MatrixXd> y0({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)}, {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});
    Y<MatrixXd> yl = y0;
    for (int i = 0; i <= match; ++i) {
        if (i > 0)
            yl = M[i - 1] * yl;
        yl = sectors[i]->propagate(E, yl, true);
    }
    Y<MatrixXd> yr = y0;
    for (int i = sectorCount - 1; i > match; --i) {
        yr = sectors[i]->propagate(E, yr, false);
        yr = M[i - 1].transpose() * yr;
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
pair<double, double> SEBase<n>::calculateError(
        double E, const std::function<bool(std::pair<double, double>, std::pair<double, double>)> &sorter) const {
    vector<pair<double, double>> *errors = sortedErrors(E, sorter);
    pair<double, double> best = (*errors)[0];
    delete errors;
    return best;
}

template<int n>
Y<MatrixXd> *SEBase<n>::computeEigenfunctionSteps(double E) const {
    int match = sectorCount / 2;

    Y<MatrixXd> *steps = new Y<MatrixXd>[sectorCount + 1];

    steps[0] = Y<MatrixXd>({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)},
                           {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});
    steps[sectorCount] = Y<MatrixXd>({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)},
                                     {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});

    for (int i = sectorCount - 1; i >= match; --i)
        steps[i] = sectors[i]->propagate(
                E, i < sectorCount - 1 ? M[i].transpose() * steps[i + 1] : steps[i + 1], false);
    Y<MatrixXd> matchRight = steps[match];

    for (int i = 0; i < match; ++i)
        steps[i + 1] = M[i] * sectors[i]->propagate(E, steps[i], true);
    Y<MatrixXd> matchLeft = steps[match];

    MatrixXd big = MatrixXd::Zero(2 * N, 2 * N);
    big << matchLeft.y[0], -matchRight.y[0], matchLeft.y[1], -matchRight.y[1];

    FullPivLU<MatrixXd> lu(big);
    MatrixXd kernel = lu.kernel();

    if (kernel.isZero(0)) {
        delete[] steps;
        return nullptr;
    }

    MatrixXd left = kernel.topRows(N);
    MatrixXd right = kernel.bottomRows(N);
    for (int i = 0; i <= match; ++i)
        steps[i] *= left;
    for (int i = sectorCount; i > match; --i)
        steps[i] *= right;

    return steps;
};


std::vector<typename dim<2>::array>
SEnD<2>::computeEigenfunction(double E, const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const {
    long nx = x.size();
    for (int i = 1; i < nx; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("SE2D::computeEigenfunction(): x has to be sorted");

    long ny = y.size();
    for (int i = 1; i < ny; ++i)
        if (y[i - 1] > y[i])
            throw runtime_error("SE2D::computeEigenfunction(): y has to be sorted");

    if (y[0] < sectors[0]->ymin || y[ny - 1] > sectors[sectorCount - 1]->ymax)
        throw runtime_error("SE2D::computeEigenfunction(): y is out of range");


    vector<ArrayXXd> result;
    Y<MatrixXd> *steps = computeEigenfunctionSteps(E);
    if (steps != nullptr) {
        int cols = (int) steps[0].y[0].cols();
        for (int i = 0; i < cols; ++i)
            result.push_back(ArrayXXd::Zero(nx, ny));

        int nextY = 0;
        int sector = 0;
        while (nextY < ny) {
            while (sector < sectorCount && y[nextY] > sectors[sector]->ymax) {
                ++sector;
                if (sector >= sectorCount)
                    throw runtime_error("SE2D::computeEigenfunction(): y is out of range");
            }

            MatrixXd B(nx, N);
            for (int j = 0; j < N; ++j)
                B.col(j) = sectors[sector]->computeEigenfunction(j, x);

            while (nextY < ny && y[nextY] <= sectors[sector]->ymax) {
                MatrixXd prod = B * sectors[sector]->propagate(E, steps[sector], y[nextY], true).y[0];
                for (int i = 0; i < cols; ++i)
                    result[i].col(nextY) = prod.col(i);
                ++nextY;
            }
        }

        delete[] steps;
    }

    return result;
}

template
class matslise::SEBase<2>;

template
class matslise::SEBase<3>;