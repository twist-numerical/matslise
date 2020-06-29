#include "../matslise.h"
#include "../util/legendre.h"
#include "../util/quadrature.h"
#include <map>

using namespace Eigen;
using namespace matslise;
using namespace std;
using namespace quadrature;

template<typename Scalar>
Array<Scalar, Dynamic, Dynamic>
eval2d(const function<Scalar(const Scalar &, const Scalar &)> &f,
       const Array<Scalar, Dynamic, 1> &x, const Array<Scalar, Dynamic, 1> &y) {
    Array<Scalar, Dynamic, Dynamic> result(x.size(), y.size());
    for (Index i = 0; i < x.size(); ++i)
        for (Index j = 0; j < x.size(); ++j)
            result(i, j) = f(x[i], y[j]);
    return result;
}


template<typename Scalar>
Matslise3D<Scalar>::Matslise3D(
        const std::function<Scalar(Scalar, Scalar, Scalar)> &potential, const matslise::Rectangle<3, Scalar> &domain,
        const SectorBuilder<Matslise3D<Scalar>> &sectorBuilder, const Scalar &tolerance)
        : AbstractMatslise3D<Scalar>(potential, domain), N(10), tolerance(tolerance) {
    grid_x = lobatto::grid<Scalar>(ArrayXs::LinSpaced(101, domain.template getMin<0>(), domain.template getMax<0>()));
    grid_y = lobatto::grid<Scalar>(ArrayXs::LinSpaced(101, domain.template getMin<1>(), domain.template getMax<1>()));

    auto sectorsBuild = sectorBuilder(this, domain.min, domain.max);
    sectors = std::move(sectorsBuild.sectors);
    matchIndex = sectorsBuild.matchIndex;
    Index sectorCount = sectors.size();

    M.reserve(sectorCount - 1);
    for (int k = 0; k < sectorCount - 1; ++k) {
        MatrixXs r(N, N);

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r(i, j) = lobatto::quadrature<Scalar>(
                        grid_x, grid_y,
                        sectors[k]->eigenfunctions_grid[j] *
                        sectors[k + 1]->eigenfunctions_grid[i]);

        cout << "M:" << endl;
        cout << r << endl;
        cout << endl;
        M.push_back(move(r));
    }
}

template<typename Scalar>
Matslise3D<Scalar>::~Matslise3D() {
    for (auto &sector : sectors)
        delete sector;
}

template<typename Scalar>
typename Matscs<Scalar>::Sector *initializeMatscs(const typename Matslise3D<Scalar>::Sector &sector) {
    const int &N = sector.matslise3d->N;
    auto &matslise3d = sector.matslise3d;
    typedef typename Matslise3D<Scalar>::MatrixXs MatrixXs;
    typedef typename Matslise3D<Scalar>::ArrayXXs ArrayXXs;
    return new typename Matscs<Scalar>::Sector(
            legendre::getCoefficients<MATSCS_N, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Scalar>(
                    [&](Scalar z) -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> {
                        MatrixXs dV(N, N);

                        ArrayXXs vDiff = eval2d<Scalar>([&](const Scalar &x, const Scalar &y) {
                            return matslise3d->potential(x, y, z);
                        }, matslise3d->grid_x, matslise3d->grid_y) - sector.vbar;

                        for (int i = 0; i < N; ++i) {
                            for (int j = 0; j <= i; ++j) {
                                dV(i, j) = lobatto::quadrature<Scalar>(
                                        matslise3d->grid_x, matslise3d->grid_y,
                                        sector.eigenfunctions_grid[i] * vDiff * sector.eigenfunctions_grid[j]);
                                if (j < i) dV(j, i) = dV(i, j);
                            }
                            dV(i, i) += sector.eigenvalues[i];
                        }

                        cout << "Z: " << z << " dV:" << endl;
                        cout << dV << "\n" << endl;
                        return dV;
                    }, sector.min, sector.max),
            sector.min, sector.max, sector.backward);
}

template<typename Scalar>
Matslise3D<Scalar>::Sector::Sector(
        const Matslise3D<Scalar> *matslise3d, const Scalar &zmin, const Scalar &zmax, bool backward)
        : matslise3d(matslise3d), min(zmin), max(zmax), backward(backward) {

    zbar = (zmax + zmin) / 2;
    function<Scalar(const Scalar &, const Scalar &)> vbar_fun = [matslise3d, this](
            const Scalar &x, const Scalar &y) -> Scalar {
        return matslise3d->potential(x, y, zbar);
    };

    matslise2d = std::make_shared<Matslise2D<Scalar>>(
            vbar_fun, matslise3d->domain.sub,
            Options2<Scalar>().tolerance(matslise3d->tolerance)
                    .nested(Options1<Scalar>().tolerance(matslise3d->tolerance)));

    cout << "Seeking: " << zbar << endl;
    eigenvalues = matslise2d->eigenvaluesByIndex(0, matslise3d->N);
    cout << "found: " << zbar << endl;
    for (auto &E : eigenvalues)
        cout << " " << E;
    cout << endl;
    if (static_cast<int>(eigenvalues.size()) != matslise3d->N) {
        throw std::runtime_error("Matlise3D: not enough basis-functions found on a sector");
    }
    eigenfunctions.reserve(matslise3d->N);
    eigenfunctions_grid.reserve(matslise3d->N);

    Index i = 0;
    for (const Scalar &E : eigenvalues) {
        for (auto &f :  matslise2d->eigenfunction(E)) {
            eigenfunctions.push_back(f);
            eigenfunctions_grid.push_back(move(f(matslise3d->grid_x, matslise3d->grid_y)));

            if (++i == matslise3d->N)
                break;
        }
        if (i == matslise3d->N)
            break;
    }
    vbar = eval2d<Scalar>([&](const Scalar &x, const Scalar &y) {
        return matslise3d->potential(x, y, zbar);
    }, matslise3d->grid_x, matslise3d->grid_y);

    matscs = initializeMatscs<Scalar>(*this);
}

template<typename Scalar>
Matslise3D<Scalar>::Sector::~Sector() {
    delete matscs;
}

template<typename Scalar>
template<int r>
Y<Scalar, Eigen::Dynamic, r>
Matslise3D<Scalar>::Sector::propagate(
        const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b, bool use_h) const {
    return matscs->propagateColumn(
            E, y0, a < min ? min : a > max ? max : a, b < min ? min : b > max ? max : b, use_h);
}

template<typename Scalar>
typename Matslise3D<Scalar>::MatrixXs Matslise3D<Scalar>::conditionY(Y<Scalar, Dynamic> &y) const {
    MatrixXs U = y.getY(0).partialPivLu().matrixLU();
    U.template triangularView<StrictlyLower>().setZero();
    U.template triangularView<Upper>().
            template solveInPlace<OnTheRight>(y.y);
    U.template triangularView<Upper>().
            template solveInPlace<OnTheRight>(y.dy);
    return U;
}

template<typename Scalar>
pair<typename Matslise3D<Scalar>::MatrixXs, typename Matslise3D<Scalar>::MatrixXs>
Matslise3D<Scalar>::matchingErrorMatrix(const Scalar &E, bool use_h) const {
    Y<Scalar, Dynamic> yl = Y<Scalar, Dynamic>::Dirichlet(N);
    for (int i = 0; i <= matchIndex; ++i) {
        yl = M[i] * sectors[i]->propagate(E, yl, sectors[i]->min, sectors[i]->max, use_h);
        conditionY(yl);
    }
    Y<Scalar, Dynamic> yr = sectors.back()->propagate(
            E, Y<Scalar, Dynamic>::Dirichlet(N), sectors.back()->max, sectors.back()->min, use_h);
    conditionY(yr);
    for (int i = sectors.size() - 2; i > matchIndex; --i) {
        yr = sectors[i]->propagate(E, (MatrixXs)(M[i].transpose()) * yr, sectors[i]->max, sectors[i]->min, use_h);
        conditionY(yr);
    }


    ColPivHouseholderQR<MatrixXs> left_solver(yl.getY(0).transpose());
    ColPivHouseholderQR<MatrixXs> right_solver(yr.getY(0).transpose());
    MatrixXs Ul = left_solver.solve(yl.getY(1).transpose()).transpose();
    MatrixXs Ur = right_solver.solve(yr.getY(1).transpose()).transpose();
    return make_pair(
            Ul - Ur,
            left_solver.solve((yl.getdY(1) - Ul * yl.getdY(0)).transpose()).transpose()
            - right_solver.solve((yr.getdY(1) - Ur * yr.getdY(0)).transpose()).transpose()
    );

}

template<typename Scalar>
vector<pair<Scalar, Scalar>> Matslise3D<Scalar>::matchingErrors(const Scalar &E, bool use_h) const {
    pair<MatrixXs, MatrixXs> error_matrix = matchingErrorMatrix(E, use_h);
    EigenSolver<MatrixXs> solver(N);

    solver.compute(error_matrix.first, true);

    multimap<Scalar, int> rightMap;
    Array<Scalar, Dynamic, 1> eigenvaluesRight = solver.eigenvalues().array().real();
    for (int i = 0; i < N; ++i)
        rightMap.insert({eigenvaluesRight[i], i});
    MatrixXs right = solver.eigenvectors().real();

    multimap<Scalar, int> leftMap;
    solver.compute(error_matrix.first.transpose(), true);
    Array<Scalar, Dynamic, 1> eigenvaluesLeft = solver.eigenvalues().array().real();
    for (int i = 0; i < N; ++i)
        leftMap.insert({eigenvaluesLeft[i], i});
    MatrixXs left = solver.eigenvectors().transpose().real();


    vector<pair<Scalar, Scalar>> errors;
    for (auto leftI = leftMap.begin(), rightI = rightMap.begin();
         leftI != leftMap.end() && rightI != rightMap.end();
         ++leftI, ++rightI) {
        int &li = leftI->second;
        int &ri = rightI->second;
        errors.push_back({eigenvaluesLeft[li],
                          (left.row(li) * error_matrix.second * right.col(ri) /
                           (left.row(li) * right.col(ri)))[0]});
    }

    return errors;
}

template<typename Scalar>
Scalar Matslise3D<Scalar>::eigenvalue(const Scalar &_E, bool use_h) const {
    const Scalar tolerance = 1e-9;
    const Scalar minTolerance = 1e-5;
    const int maxIterations = 30;

    Scalar E = _E;
    Scalar error, derror;
    int i = 0;
    do {
        tie(error, derror) = matchingError(E, use_h);
        E -= error / derror;
        ++i;
    } while (i < maxIterations && abs(error) > tolerance);

    if (abs(error) > minTolerance)
        return NAN;
    return E;
}


template<typename Scalar>
Scalar Matslise3D<Scalar>::eigenvalueError(const Scalar &E) const {
    return abs(E - eigenvalue(E, false));
}

template<typename Scalar>
pair<Scalar, Scalar>
Matslise3D<Scalar>::matchingError(const Scalar &E, bool use_h) const {
    vector<pair<Scalar, Scalar>> errors = matchingErrors(E, use_h);
    return *min_element(errors.begin(), errors.end(), [](
            const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) -> bool {
        if (abs(a.first) > 100 || abs(b.first) > 100)
            return abs(a.first) < abs(b.first);
        return abs(a.first / a.second) < abs(b.first / b.second);
    });
}

#define INSTANTIATE_PROPAGATE(Scalar, r) \
template Y<Scalar, Dynamic, r> \
Matslise3D<Scalar>::Sector::propagate<r>( \
        const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b, bool use_h) const;

template<typename Scalar>
Scalar Matslise3D<Scalar>::Sector::error() const {
    return matscs->error();
}


template<typename Scalar>
typename Matslise3D<Scalar>::Sector *Matslise3D<Scalar>::Sector::refine(
        const Matslise3D<Scalar> *problem, const Scalar &_min, const Scalar &_max, bool _backward) const {
    Scalar h = _max - _min;
    if (backward != _backward || zbar < _min + h / 3 || zbar > _max - h / 3) {
        return new Sector(problem, _min, _max, _backward);
    }

    // std::cout << "refining sector 2d" << std::endl;
    auto sector = new Sector(problem);
    sector->min = _min;
    sector->max = _max;
    sector->backward = _backward;
    sector->zbar = zbar;
    sector->matslise2d = matslise2d;
    sector->eigenfunctions = eigenfunctions;
    sector->eigenvalues = eigenvalues;
    sector->eigenfunctions_grid = eigenfunctions_grid;
    sector->matscs = initializeMatscs<Scalar>(*sector);
    return sector;
}

#define INSTANTIATE_MORE(Scalar) \
INSTANTIATE_PROPAGATE(Scalar, 1) \
INSTANTIATE_PROPAGATE(Scalar, -1)

#include "../util/instantiate.h"