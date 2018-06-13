//
// Created by toon on 6/13/18.
//

#include <iostream>
#include "../se2d.h"
#include "../lobatto.h"

#define N 12

using namespace se2d;
using namespace std;

#define GRID_POINTS 40

template<typename FROM, typename TO>
Array<TO, Dynamic, 1> apply(const Array<FROM, Dynamic, 1> &x, function<TO(FROM)> f) {
    Array<TO, Dynamic, 1> r(x.size());
    for (int i = 0; i < r.size(); ++i)
        r[i] = f(x[i]);
    return r;
}

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

    cout << M << endl<< endl;

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

matscs::Y SE2D::propagate(double E, double y, bool forward) const {
    matscs::Y c(MatrixXd::Zero(N, N), MatrixXd::Identity(N, N));
    if (forward) {
        int i = 0;

        Sector *sector = sectors[i];
        while (i < sectorCount && y > sector->ymax) {
            if (i > 0) {
                c.y = M[i - 1] * c.y;
                c.dy = M[i - 1] * c.dy;
            }
            c = sector->propagate(E, c, true);
            sector = sectors[++i];
        }

        if (i < sectorCount && y > sector->ymin) {
            if (i > 0) {
                c.y = M[i - 1] * c.y;
                c.dy = M[i - 1] * c.dy;
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


Sector::Sector(SE2D *se2d, double ymin, double ymax, int sectorCount) : se2d(se2d), ymin(ymin), ymax(ymax) {
    static const matslise::Y y0 = matslise::Y({0, 1});

    const double ybar = (ymax + ymin) / 2;
    function<double(double)> vbar_fun = [se2d, ybar](double x) -> double { return se2d->V(x, ybar); };
    vbar = apply(se2d->xGrid, vbar_fun);
    matslise = new Matslise(vbar_fun, se2d->xmin, se2d->xmax, sectorCount);

    vector<tuple<unsigned int, double>> *index_eigv = matslise->computeEigenvaluesByIndex(0, N, y0, y0);
    eigenvalues = new double[N];
    eigenfunctions = new ArrayXd[N];
    for (int i = 0; i < N; ++i) {
        double E = get<1>((*index_eigv)[i]);
        eigenvalues[i] = E;
        Array<matslise::Y, Dynamic, 1> func = matslise->computeEigenfunction(E, y0, y0, se2d->xGrid);
        eigenfunctions[i] = ArrayXd(func.size());
        for (int j = 0; j < func.size(); ++j)
            eigenfunctions[i][j] = func[j].y[0];
        eigenfunctions[i] /= sqrt(lobatto::quadrature(se2d->xGrid, eigenfunctions[i] * eigenfunctions[i]));
    }

    matscs = new Matscs([this](double y) -> MatrixXd { return this->calculateDeltaV(y); }, N, ymin, ymax, 3);

    delete index_eigv;
}

Sector::~Sector() {
    delete matslise;
    delete[] eigenvalues;
    delete[] eigenfunctions;
}

matscs::Y Sector::propagate(double E, const matscs::Y &c, bool forward) const {
    return propagate(E, c, forward ? ymax : ymin, forward);
}

matscs::Y Sector::propagate(double E, const matscs::Y &c, double y, bool forward) const {
    return matscs->propagate(E, c, forward ? ymin : ymax, y);
}

Eigen::MatrixXd Sector::calculateDeltaV(double y) const {
    Eigen::MatrixXd dV(N, N);

    ArrayXd vDiff =
            apply<double, double>(se2d->xGrid, [this, y](double x) -> double { return this->se2d->V(x, y); }) - vbar;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            dV(i, j) = lobatto::quadrature(se2d->xGrid, eigenfunctions[i] * eigenfunctions[j] * vDiff);
            if(j < i) dV(j, i) = dV(i, j);
        }
        dV(i, i) += eigenvalues[i];
    }

    return dV;
}
