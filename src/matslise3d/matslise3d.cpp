#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/scoped_timer.h"
#include <map>

using namespace Eigen;
using namespace matslise;
using namespace std;
using namespace quadrature;

template<typename Scalar>
Matslise3D<Scalar>::Matslise3D(
        const std::function<Scalar(Scalar, Scalar, Scalar)> &potential,
        const matslise::Rectangle<Scalar, 3> &domain, const Config &config)
        : AbstractMatslise3D<Scalar>(potential, domain), config(config) {
    MATSLISE_SCOPED_TIMER("3D constructor");
    grid_x = lobatto::grid<Scalar>(ArrayXs::LinSpaced(101, domain.template min<0>(), domain.template max<0>()));
    grid_y = lobatto::grid<Scalar>(ArrayXs::LinSpaced(101, domain.template min<1>(), domain.template max<1>()));

    auto sectorsBuild = sector_builder::getOrAutomatic<Matslise3D<Scalar>, true>(config.zSectorBuilder, config.tolerance)
            (this, domain.template min<2>(), domain.template max<2>());
    sectors = std::move(sectorsBuild.sectors);
    matchIndex = sectorsBuild.matchIndex;
    Index sectorCount = sectors.size();
    cout << "sectors build" << endl;

    M.reserve(sectorCount - 1);
    const Index &N = config.xyBasisSize;
    for (int k = 0; k < sectorCount - 1; ++k) {
        MatrixXs &r = M.emplace_back(N, N);

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r(i, j) = lobatto::quadrature<Scalar>(
                        grid_x, grid_y,
                        sectors[k]->eigenfunctions_grid[j] *
                        sectors[k + 1]->eigenfunctions_grid[i]);
        cout << "\n\nM " << k << endl;
        cout << r * r.transpose() - MatrixXs::Identity(N, N) << endl;
    }
}

template<typename Scalar>
Matslise3D<Scalar>::~Matslise3D() {
    for (auto &sector : sectors)
        delete sector;
}

template<typename Scalar>
typename Matslise3D<Scalar>::MatrixXs Matslise3D<Scalar>::conditionY(Y<Scalar, Dynamic> &y) const {
    MATSLISE_SCOPED_TIMER("3D conditionY");
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
    MATSLISE_SCOPED_TIMER("3D matchingErrorMatrix");
    const Index &N = config.xyBasisSize;
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
    MATSLISE_SCOPED_TIMER("3D matchingErrors");
    pair<MatrixXs, MatrixXs> error_matrix = matchingErrorMatrix(E, use_h);
    const Index &N = config.xyBasisSize;
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
    MATSLISE_SCOPED_TIMER("3D eigenvalue");
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

#include "../util/instantiate.h"