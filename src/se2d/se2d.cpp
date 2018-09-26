//
// Created by toon on 6/13/18.
//

#include <iostream>
#include "../se2d.h"
#include "../util/lobatto.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::se2d_util;
using namespace std;


SE2D::SE2D(function<double(double, double)> V,
           double xmin, double xmax, double ymin, double ymax, int xSectorCount, int ySectorCount) :
        V(V), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), sectorCount(ySectorCount) {
    sectors = new Sector *[sectorCount];


    ArrayXd xs(GRID_POINTS);
    for (int i = 0; i < GRID_POINTS; ++i)
        xs[i] = xmin + (xmax - xmin) * i / (GRID_POINTS - 1);
    xGrid = lobatto::grid(xs);

    double h = (ymax - ymin) / sectorCount;
    for (int i = 0; i < sectorCount; ++i)
        sectors[i] = new Sector(this, ymin + i * h, ymin + (i + 1) * h, xSectorCount);

    M = new MatrixXd[sectorCount - 1];
    for (int i = 0; i < sectorCount - 1; ++i)
        M[i] = calculateM(i);
}

MatrixXd SE2D::calculateM(int k) const {
    MatrixXd M(N, N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            M(i, j) = lobatto::quadrature(xGrid, sectors[k]->eigenfunctions[j] * sectors[k + 1]->eigenfunctions[i]);

    return M;
}


SE2D::~SE2D() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
    delete[] M;
}

MatrixXd SE2D::calculateErrorMatrix(double E) const {
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
    return yl.y.y * yl.y.x.inverse() - yr.y.y * yr.y.x.inverse();
}

double SE2D::calculateError(double E) const {
    ArrayXcd eigenvalues = calculateErrorMatrix(E).eigenvalues().array();
    ArrayXcd::Index index;
    eigenvalues.abs2().minCoeff(&index);
    return eigenvalues[index].real();
}

Y<MatrixXd> *SE2D::computeEigenfunctionSteps(double E) const {
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
    big << matchLeft.y.x, -matchRight.y.x, matchLeft.y.y, -matchRight.y.y;

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

vector<ArrayXXd>
SE2D::computeEigenfunction(double E, const ArrayXd &x, const ArrayXd &y) const {
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
        int cols = steps[0].y.x.cols();
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
                MatrixXd prod = B * sectors[sector]->propagate(E, steps[sector], y[nextY], true).y.x;
                for (int i = 0; i < cols; ++i)
                    result[i].col(nextY) = prod.col(i);
                ++nextY;
            }
        }

        delete[] steps;
    }

    return result;
}


const Sector &SE2D::getSector(double y) const {
    for (int i = 0; i < sectorCount; ++i)
        if (y < sectors[i]->ymax)
            return *sectors[i];
    throw runtime_error("SE2D::getSector(): no sector found");
}