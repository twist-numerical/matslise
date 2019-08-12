#include "../matslise.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;

#define cec_cce(y) ((y).getdY(0).transpose()*(y).getY(1) - (y).getY(0).transpose()*(y).getdY(1))

template<typename Scalar>
vector<Y<Scalar, Dynamic>> SE2D<Scalar>::computeEigenfunctionSteps(const Scalar &E) const {
    Y<Scalar, Dynamic> *steps = new Y<Scalar, Dynamic>[sectorCount + 1];

    steps[0] = Y<Scalar, Dynamic>::Dirichlet(N);
    steps[sectorCount] = Y<Scalar, Dynamic>::Dirichlet(N);

    MatrixXs normRight = MatrixXs::Zero(N, N);
    int matchIndex = 0;
    for (int i = sectorCount - 1; sectors[i]->min > match; --i) {
        Y<Scalar, Dynamic> next
                = i < sectorCount - 1 ? (MatrixXs)(M[i].transpose()) * steps[i + 1] : steps[i + 1];
        steps[i] = sectors[i]->propagate(E, next, sectors[i]->max, sectors[i]->min, true);
        normRight += cec_cce(next) - cec_cce(steps[i]);
        matchIndex = i;
    }
    Y<Scalar, Dynamic> matchRight = steps[matchIndex];

    MatrixXs normLeft = MatrixXs::Zero(N, N);
    for (int i = 0; i < matchIndex; ++i) {
        Y<Scalar, Dynamic> next = sectors[i]->propagate(E, steps[i], sectors[i]->min, sectors[i]->max, true);
        steps[i + 1] = M[i] * next;
        normLeft += cec_cce(next) - cec_cce(steps[i]);
    }
    Y<Scalar, Dynamic> matchLeft = steps[matchIndex];

    ColPivHouseholderQR<MatrixXs> left_solver(matchLeft.getY(0).transpose());
    ColPivHouseholderQR<MatrixXs> right_solver(matchRight.getY(0).transpose());
    MatrixXs Ul = left_solver.solve(matchLeft.getY(1).transpose()).transpose();
    MatrixXs Ur = right_solver.solve(matchRight.getY(1).transpose()).transpose();

    FullPivLU<MatrixXs> lu(Ul - Ur);
    lu.setThreshold(1e-4);

    vector<Y<Scalar, Dynamic>> elements;
    if (lu.dimensionOfKernel() > 0) {
        elements.resize(static_cast<size_t>(sectorCount + 1));
        MatrixXs kernel = lu.kernel();

        MatrixXs left = matchLeft.getY(0).colPivHouseholderQr().solve(kernel);
        MatrixXs right = matchRight.getY(0).colPivHouseholderQr().solve(kernel);
        MatrixXs scaling = (left.transpose() * normLeft * left
                            + right.transpose() * normRight * right).diagonal();
        scaling = scaling.unaryExpr([](Scalar s) { return s < 0 ? 1 : 1. / sqrt(s); });
        left *= scaling.asDiagonal();
        right *= scaling.asDiagonal();
        for (int i = 0; i <= matchIndex; ++i)
            elements[static_cast<size_t>(i)] = steps[i] * left;
        for (int i = sectorCount; i > matchIndex; --i)
            elements[static_cast<size_t>(i)] = steps[i] * right;
    }
    delete[] steps;
    return elements;
}

template<typename Scalar>
std::vector<typename SE2D<Scalar>::ArrayXXs>
SE2D<Scalar>::computeEigenfunction(
        const Scalar &E, const typename SE2D<Scalar>::ArrayXs &x, const typename SE2D<Scalar>::ArrayXs &y) const {

    Eigen::Index nx = x.size();
    for (Eigen::Index i = 1; i < nx; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("SE2D::computeEigenfunction(): x has to be sorted");

    Eigen::Index ny = y.size();
    for (Eigen::Index i = 1; i < ny; ++i)
        if (y[i - 1] > y[i])
            throw runtime_error("SE2D::computeEigenfunction(): y has to be sorted");

    if (y[0] < sectors[0]->min || y[ny - 1] > sectors[sectorCount - 1]->max)
        throw runtime_error("SE2D::computeEigenfunction(): y is out of range");


    vector<ArrayXXs> result;
    vector<Y<Scalar, Dynamic>> steps = computeEigenfunctionSteps(E);
    if (steps.size() > 0) {
        Eigen::Index cols = steps[0].getY(0).cols();
        for (Eigen::Index i = 0; i < cols; ++i)
            result.push_back(ArrayXXs::Zero(nx, ny));

        Eigen::Index nextY = 0;
        int sector = 0;
        while (nextY < ny) {
            while (sector < sectorCount && y[nextY] > sectors[sector]->max) {
                ++sector;
                if (sector >= sectorCount)
                    throw runtime_error("SE2D::computeEigenfunction(): y is out of range");
            }

            MatrixXs B(nx, N);
            for (int j = 0; j < N; ++j)
                B.col(j) = sectors[sector]->computeEigenfunction(j, x);

            while (nextY < ny && y[nextY] <= sectors[sector]->max) {
                MatrixXs prod = B * sectors[sector]->propagate(
                        E, steps[static_cast<size_t>(sector)], sectors[sector]->min, y[nextY], true).getY(0);
                for (int i = 0; i < cols; ++i)
                    result[static_cast<unsigned>(i)].col(nextY) = prod.col(i);
                ++nextY;
            }
        }
    }

    return result;
}

template<typename Scalar>
vector<function<Scalar(Scalar, Scalar)>> SE2D<Scalar>::eigenfunctionCalculator(const Scalar &E) const {
    shared_ptr<vector<Y<Scalar, Dynamic>>> steps
            = make_shared<vector<Y<Scalar, Dynamic>>>(move(computeEigenfunctionSteps(E)));
    vector<function<Scalar(Scalar, Scalar)>> result;
    auto bases = make_shared<vector<function<ArrayXs(Scalar)>>>(sectorCount);
    for (int i = 0; i < sectorCount; ++i)
        (*bases)[static_cast<size_t>(i)] = sectors[i]->basisCalculator();

    if (steps->size() > 0) {
        int cols = (int) steps->at(0).getY(0).cols();
        for (int column = 0; column < cols; ++column) {
            result.push_back([this, steps, E, column, bases](Scalar x, Scalar y) -> Scalar {
                unsigned long sectorIndex;
                {
                    unsigned long a = 0;
                    unsigned long b = (unsigned long) this->sectorCount;
                    while (a + 1 < b) {
                        unsigned long c = (a + b) / 2;
                        if (y < this->sectors[c]->min)
                            b = c;
                        else
                            a = c;
                    }
                    sectorIndex = a;
                }
                const SE2D<Scalar>::Sector *sector = this->sectors[sectorIndex];

                return sector->propagate(
                        E, (*steps)[sectorIndex].col(column), sector->min, y, true).getY(0).dot(
                        (*bases)[sectorIndex](x).matrix()
                );
            });
        }
    }
    return result;
}

#include "../util/instantiate.h"