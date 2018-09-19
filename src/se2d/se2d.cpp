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

MatrixXd SE2D::calculateError(double E) const {
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


Y<MatrixXd> SE2D::propagate(double E, double y, bool forward) const {
    Y<MatrixXd> c({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)}, {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});
    if (forward) {
        int i = 0;

        Sector *sector = sectors[i];
        while (i < sectorCount && y > sector->ymax) {
            if (i > 0) {
                c.y = M[i - 1] * c.y;
            }
            c = sector->propagate(E, c, true);
            sector = sectors[++i];
        }

        if (i < sectorCount && y > sector->ymin) {
            if (i > 0) {
                c.y = M[i - 1] * c.y;
            }
            c = sector->propagate(E, c, y, true);
        }
    } else {
        int i = sectorCount - 1;

        Sector *sector = sectors[i];
        while (i >= 0 && y < sector->ymin) {
            c = sector->propagate(E, c, false);
            sector = sectors[--i];
        }

        if (i >= 0 && y < sector->ymax)
            c = sector->propagate(E, c, y, false);
    }

    return c;
}


ArrayXXd
SE2D::computeEigenfunction(double E, const ArrayXd &x, const ArrayXd &y) const {
    long nx = x.size();
    for (int i = 1; i < nx; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("SE2D::computeEigenfunction(): x has to be sorted");

    long ny = y.size();
    for (int i = 1; i < ny; ++i)
        if (y[i - 1] > y[i])
            throw runtime_error("SE2D::computeEigenfunction(): y has to be sorted");

    int match = sectorCount / 2;

    ArrayXXd result(nx, ny);
    Array<Y<MatrixXd>, Dynamic, 1> steps(sectorCount + 1);

    steps(0) = Y<MatrixXd>({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)}, {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});
    steps(sectorCount) = Y<MatrixXd>({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)}, {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});

    for (int i = 1; i <= match; ++i)
        steps(i) = sectors[i]->propagate(E, steps(i - 1), true);

    for (int i = sectorCount - 1; i > match; --i)
        steps(i) = sectors[i]->propagate(E, steps(i + 1), false);

    MatrixXd big = MatrixXd::Zero(2 * N, 2 * N);
    big << steps(match).y.x, -steps(match + 1).y.x, steps(match).y.y, -steps(match + 1).y.y;

    FullPivLU<MatrixXd> lu(big);
    MatrixXd A_null_space = lu.kernel();
    cout << A_null_space << endl;

    return result;
}