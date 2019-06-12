//
// Created by toon on 6/13/18.
//

#include <iostream>
#include <map>
#include "../matslise.h"
#include "../util/lobatto.h"
#include "../util/find_sector.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;


template<typename Scalar>
Array<Scalar, Dynamic, 1> getGrid(const Scalar &min, const Scalar &max, int count) {
    Array<Scalar, Dynamic, 1> points(count);
    for (int i = 0; i < count; ++i)
        points[i] = min + (max - min) * i / (count - 1);
    return points;
}

template<typename Scalar>
SE2D<Scalar>::SE2D(const function<Scalar(const Scalar &, const Scalar &)> &V, const Rectangle<2, Scalar> &domain,
                   const Options2<Scalar> &_options) :
        V(V), domain(domain), N(_options._N),
        options(_options) {
    grid = lobatto::grid<Scalar>(getGrid(domain.getMin(0), domain.getMax(0), options._gridPoints));
    options._builder->build(this, domain.min, domain.max);

    M = new MatrixXs[sectorCount - 1];
    for (int i = 0; i < sectorCount - 1; ++i)
        M[i] = calculateM(i);
}

template<typename Scalar>
typename SE2D<Scalar>::MatrixXs SE2D<Scalar>::calculateM(int k) const {
    MatrixXs M(N, N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            M(i, j) = lobatto::quadrature<Scalar>(grid,
                                                  sectors[k]->eigenfunctions[j] * sectors[k + 1]->eigenfunctions[i]);

    return M;
}

template<typename Scalar>
SE2D<Scalar>::~SE2D() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
    delete[] M;
}

template<typename Scalar>
pair<typename SE2D<Scalar>::MatrixXs, typename SE2D<Scalar>::MatrixXs>
SE2D<Scalar>::calculateErrorMatrix(const Scalar &E) const {
    Y<Scalar, Dynamic> y0 = Y<Scalar, Dynamic>::Dirichlet(N);
    Y<Scalar, Dynamic> yl = y0;
    for (int i = 0; sectors[i]->max <= match; ++i) {
        yl = M[i] * sectors[i]->propagate(E, yl, sectors[i]->min, sectors[i]->max, true);
        //yl *= (yl.getY(0).colwise().norm()).cwiseInverse().asDiagonal();
    }
    Y<Scalar, Dynamic> yr = sectors[sectorCount - 1]->propagate(
            E, y0, sectors[sectorCount - 1]->max, sectors[sectorCount - 1]->min, true);
    for (int i = sectorCount - 2; sectors[i]->min >= match; --i) {
        yr = sectors[i]->propagate(E, (MatrixXs)(M[i].transpose()) * yr, sectors[i]->max, sectors[i]->min, true);
        //yr *= (yr.getY(0).colwise().norm()).cwiseInverse().asDiagonal();
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
Y<Scalar, Dynamic>
SE2D<Scalar>::propagate(const Scalar &E, const Y<Scalar, Dynamic> &y0, const Scalar &a, const Scalar &b,
                        bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Scalar, Dynamic> y = y0;
    int sectorIndex = find_sector<SE2D<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        if (direction == -1 && sectorIndex < sectorCount - 1)
            y = (MatrixXs)(M[sectorIndex].transpose()) * y;
        y = sector->propagate(E, y, a, b, use_h);
        if (direction == 1)
            y = M[sectorIndex] * y;
        sectorIndex += direction;
    } while (!sector->contains(b));
    return y;
}

template<typename Scalar>
vector<pair<Scalar, Scalar>> *SE2D<Scalar>::calculateErrors(const Scalar &E) const {
    pair<MatrixXs, MatrixXs> error_matrix = calculateErrorMatrix(E);
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


    auto errors = new vector<pair<Scalar, Scalar>>();
    for (auto leftI = leftMap.begin(), rightI = rightMap.begin();
         leftI != leftMap.end() && rightI != rightMap.end();
         ++leftI, ++rightI) {
        int &li = leftI->second;
        int &ri = rightI->second;
        errors->push_back({eigenvaluesLeft[li],
                           (left.row(li) * error_matrix.second * right.col(ri) /
                            (left.row(li) * right.col(ri)))[0]});
    }

    return errors;
}

template<typename Scalar>
vector<pair<Scalar, Scalar>> *
SE2D<Scalar>::sortedErrors(
        const Scalar &E,
        const std::function<bool(pair<Scalar, Scalar>, pair<Scalar, Scalar>)> &sorter) const {
    vector<pair<Scalar, Scalar>> *errors = calculateErrors(E);

    sort(errors->begin(), errors->end(), sorter);

    return errors;
}

template<typename Scalar>
Scalar SE2D<Scalar>::findEigenvalue(
        const Scalar &_E, const Scalar &tolerance, int maxIterations, const Scalar &minTolerance) const {
    Scalar E = _E;
    Scalar error, derror;
    int i = 0;
    do {
        tie(error, derror) = calculateError(E, &NEWTON_RAPHSON_SORTER<Scalar>);
        E -= error / derror;
        ++i;
    } while (i < maxIterations && abs(error) > tolerance);

    if (abs(error) > minTolerance)
        return NAN;
    return E;
}

template<typename Scalar>
pair<Scalar, Scalar> SE2D<Scalar>::calculateError(
        const Scalar &E,
        const std::function<bool(std::pair<Scalar, Scalar>, std::pair<Scalar, Scalar>)> &sorter) const {
    vector<pair<Scalar, Scalar>> *errors = sortedErrors(E, sorter);
    pair<Scalar, Scalar> best = (*errors)[0];
    delete errors;
    return best;
}

#include "../util/instantiate.h"