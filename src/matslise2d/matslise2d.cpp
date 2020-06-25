#include <iostream>
#include <map>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/find_sector.h"

using namespace Eigen;
using namespace matslise;
using namespace std;
using namespace quadrature;


template<typename Scalar>
Array<Scalar, Dynamic, 1> getGrid(const Scalar &min, const Scalar &max, int count) {
    Array<Scalar, Dynamic, 1> points(count);
    for (int i = 0; i < count; ++i)
        points[i] = min + (max - min) * i / (count - 1);
    return points;
}

template<typename Scalar>
Matslise2D<Scalar>::Matslise2D(const function<Scalar(Scalar, Scalar)> &potential,
                               const matslise::Rectangle<2, Scalar> &domain,
                               const Options2<Scalar> &_options):
        AbstractMatslise2D<Scalar>(potential, domain), N(_options._N), options(_options) {
    dirichletBoundary = Y<Scalar, Eigen::Dynamic>::Dirichlet(N);
    auto sectorsBuild = options._builder(this, domain.min, domain.max);
    sectors = std::move(sectorsBuild.sectors);
    matchIndex = sectorsBuild.matchIndex;
    sectorCount = sectors.size();

/*
    vector<MatrixXs> altM;
    ArrayXs grid = lobatto::grid<Scalar>(getGrid(domain.getMin(0), domain.getMax(0), options._gridPoints));
    altM.reserve(sectorCount - 1);
    ArrayXXs basis[2]{ArrayXXs(), sectors[0]->template basis<false>(grid)};
    for (int k = 0; k < sectorCount - 1; ++k) {
        basis[k & 1] = sectors[k + 1]->template basis<false>(grid);

        altM.emplace_back(N, N);

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                altM.back()(i, j) = lobatto::quadrature<Scalar>(
                        grid, basis[k & 1].col(i) * basis[1 - (k & 1)].col(j));
    }
*/

    M.reserve(sectorCount - 1);

    map<pair<int, int>, ArrayXs> prev;
    map<pair<int, int>, ArrayXs> next;
    for (int k = 0; k < sectorCount - 1; ++k) {
        next.clear();

        M.push_back(move(
                gauss_konrod::adaptive<Scalar, MatrixXs, true>([&, k](const ArrayXs &x) {
                    ArrayXXs prevBasis = sectors[k]->template basis<false>(x);
                    ArrayXXs nextBasis = sectors[k + 1]->template basis<false>(x);

                    Array<MatrixXs, Dynamic, 1> result(x.size());
                    for (Index i = 0; i < x.size(); ++i) {
                        result(i) = nextBasis.row(i).matrix().transpose() * prevBasis.row(i).matrix();
                    }
                    return result;
                }, domain.getMin(0), domain.getMax(0), 1e-8, [&](const MatrixXs &v) {
                    return v.maxCoeff();
                })
        ));
    }
}

template<typename Scalar>
Matslise2D<Scalar>::~Matslise2D() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
}

template<typename Scalar>
typename Matslise2D<Scalar>::MatrixXs Matslise2D<Scalar>::conditionY(Y<Scalar, Dynamic> &y) const {
    MatrixXs U = y.getY(0).partialPivLu().matrixLU();
    U.template triangularView<StrictlyLower>().setZero();
    U.template triangularView<Upper>().
            template solveInPlace<OnTheRight>(y.y);
    U.template triangularView<Upper>().
            template solveInPlace<OnTheRight>(y.dy);
    return U;
}

template<typename Scalar>
pair<typename Matslise2D<Scalar>::MatrixXs, typename Matslise2D<Scalar>::MatrixXs>
Matslise2D<Scalar>::matchingErrorMatrix(const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    Y<Scalar, Dynamic> yl = yLeft;
    for (int i = 0; i <= matchIndex; ++i) {
        yl = M[i] * sectors[i]->propagate(E, yl, sectors[i]->min, sectors[i]->max, use_h);
        conditionY(yl);
    }
    Y<Scalar, Dynamic> yr = sectors[sectorCount - 1]->propagate(
            E, dirichletBoundary, sectors[sectorCount - 1]->max, sectors[sectorCount - 1]->min, use_h);
    conditionY(yr);
    for (int i = sectorCount - 2; i > matchIndex; --i) {
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
Y<Scalar, Dynamic>
Matslise2D<Scalar>::propagate(
        const Scalar &E, const Y<Scalar, Dynamic> &y0, const Scalar &a, const Scalar &b, bool use_h) const {
    if (!domain.contains(1, a) || !domain.contains(1, b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Scalar, Dynamic> y = y0;
    int sectorIndex = find_sector<Matslise2D<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        if (direction == -1 && sectorIndex < sectorCount - 1)
            y = (MatrixXs)(M[sectorIndex].transpose()) * y;
        y = sector->propagate(E, y, a, b, use_h);
        conditionY(y);
        if (direction == 1)
            y = M[sectorIndex] * y;
        sectorIndex += direction;
    } while (!sector->contains(b));
    return y;
}

template<typename Scalar>
vector<pair<Scalar, Scalar>> Matslise2D<Scalar>::matchingErrors(
        const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    pair<MatrixXs, MatrixXs> error_matrix = matchingErrorMatrix(yLeft, E, use_h);
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
pair<Scalar, Scalar>
Matslise2D<Scalar>::matchingError(const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    vector<pair<Scalar, Scalar>> errors = matchingErrors(yLeft, E, use_h);
    return *min_element(errors.begin(), errors.end(), [](
            const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) -> bool {
        if (abs(a.first) > 100 || abs(b.first) > 100)
            return abs(a.first) < abs(b.first);
        return abs(a.first / a.second) < abs(b.first / b.second);
    });
}

#include "../util/instantiate.h"