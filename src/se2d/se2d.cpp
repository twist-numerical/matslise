//
// Created by toon on 6/13/18.
//

#include <iostream>
#include "../se2d.h"
#include "../util/lobatto.h"

#define N 12

using namespace Eigen;
using namespace matslise;
using namespace matslise::se2d_sector;
using namespace std;

#define GRID_POINTS 60

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

std::tuple<double, double> SE2D::calculateError(double E) const {

    return tuple<double, double>();
}

Y<MatrixXd> SE2D::propagate(double E, double y, bool forward) const {
    Y<MatrixXd> c({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)}, {MatrixXd::Zero(N,N), MatrixXd::Zero(N,N)});
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
