//
// Created by toon on 6/13/18.
//

#include "../se2d.h"
#include "../lobatto.h"

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
}

SE2D::~SE2D() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
}


Sector::Sector(SE2D *se2d, double ymin, double ymax, int sectorCount) : se2d(se2d), ymin(ymin), ymax(ymax) {
    static const matslise::Y y0 = matslise::Y({0, 1});
    static const int N = 12;

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
        eigenfunctions[i] /= lobatto::quadrature(se2d->xGrid, eigenfunctions[i] * eigenfunctions[i]);
    }
    delete index_eigv;
}

Sector::~Sector() {
    delete matslise;
    delete[] eigenvalues;
    delete[] eigenfunctions;
}
