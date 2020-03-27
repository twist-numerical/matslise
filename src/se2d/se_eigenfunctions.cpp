#include "../matslise.h"

using namespace Eigen;
using namespace matslise;
using namespace std;

template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic> cec_cce(const Y<Scalar, Dynamic, Dynamic> &y) {
    return ((y).getdY(0).transpose() * (y).getY(1) - (y).getY(0).transpose() * (y).getdY(1));
}

template<typename Scalar>
vector<Y<Scalar, Dynamic>> SE2D<Scalar>::eigenfunctionSteps(const Scalar &E) const {
    auto *steps = new Y<Scalar, Dynamic>[sectorCount + 1];

    steps[0] = y0Left;
    steps[sectorCount] = y0Right;
    auto *U = new MatrixXs[sectorCount + 1];

    int matchIndex = 0;
    for (int i = sectorCount - 1; i >= 0 && sectors[i]->min >= match; --i) {
        const Y<Scalar, Dynamic> next
                = i < sectorCount - 1 ? (MatrixXs)(M[i].transpose()) * steps[i + 1] : steps[i + 1];
        steps[i] = sectors[i]->propagate(E, next, sectors[i]->max, sectors[i]->min, true);
        U[i + 1] = conditionY(steps[i]);
        matchIndex = i;
    }
    const Y<Scalar, Dynamic> matchRight = steps[matchIndex];

    U[0] = MatrixXs::Identity(N, N);
    for (int i = 0; i < matchIndex; ++i) {
        steps[i + 1] = M[i] * sectors[i]->propagate(E, steps[i], sectors[i]->min, sectors[i]->max, true);
        U[i + 1] = conditionY(steps[i + 1]);
    }
    const Y<Scalar, Dynamic> matchLeft = steps[matchIndex];

    ColPivHouseholderQR<MatrixXs> left_solver(matchLeft.getY(0).transpose());
    ColPivHouseholderQR<MatrixXs> right_solver(matchRight.getY(0).transpose());
    FullPivLU<MatrixXs> lu(
            left_solver.solve(matchLeft.getY(1).transpose()).transpose()
            - right_solver.solve(matchRight.getY(1).transpose()).transpose());
    lu.setThreshold(1e-4);

    vector<Y<Scalar, Dynamic>> elements;
    if (lu.dimensionOfKernel() > 0) {
        elements.resize(sectorCount + 1);
        MatrixXs kernel = lu.kernel();

        MatrixXs left = matchLeft.getY(0).colPivHouseholderQr().solve(kernel);
        MatrixXs right = matchRight.getY(0).colPivHouseholderQr().solve(kernel);
        /*
         MatrixXs scaling = (left.transpose() * normLeft * left
                            + right.transpose() * normRight * right).diagonal();
        scaling = scaling.unaryExpr([](Scalar s) { return s < 0 ? 1 : 1. / sqrt(s); });
        left *= scaling.asDiagonal();
        right *= scaling.asDiagonal();
        */

        VectorXs normalizer = (cec_cce<>(matchLeft * left) - cec_cce<>(matchRight * right)).diagonal()
                .unaryExpr([](Scalar s) { return s <= 0 ? Scalar(1) : Scalar(1) / sqrt(s); });

        for (int i = matchIndex; i >= 0; --i) {
            elements[static_cast<size_t>(i)] = steps[i] * left;
            if (i > 0) {
                U[i].template triangularView<Upper>().
                        template solveInPlace<OnTheLeft>(left);
            }
        }

        Y<Scalar, Dynamic> elementMatchRight = matchRight * right;
        for (int i = matchIndex + 1; i <= sectorCount; ++i) {
            U[i].template triangularView<Upper>().
                    template solveInPlace<OnTheLeft>(right);
            elements[static_cast<size_t>(i)] = steps[i] * right;
        }

        for (int i = 0; i <= sectorCount; ++i) {
            elements[static_cast<size_t>(i)] *= normalizer.asDiagonal();
        }
    }
    delete[] steps;
    delete[] U;
    return elements;
}

template<typename Scalar>
std::vector<typename SE2D<Scalar>::ArrayXXs>
SE2D<Scalar>::eigenfunction(
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
    vector<Y<Scalar, Dynamic>> steps = eigenfunctionSteps(E);
    if (!steps.empty()) {
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
                B.col(j) = sectors[sector]->eigenfunction(j, x);

            while (nextY < ny && y[nextY] <= sectors[sector]->max) {
                MatrixXs prod = B * sectors[sector]->propagate(
                        E, steps[static_cast<size_t>(sector)], sectors[sector]->min, y[nextY], true).getY(0);
                for (Eigen::Index i = 0; i < cols; ++i)
                    result[i].col(nextY) = prod.col(i);
                ++nextY;
            }
        }
    }

    return result;
}

template<typename Scalar>
vector<function<Scalar(Scalar, Scalar)>> SE2D<Scalar>::eigenfunctionCalculator(const Scalar &E) const {
    shared_ptr<vector<Y<Scalar, Dynamic>>> steps
            = make_shared<vector<Y<Scalar, Dynamic>>>(move(eigenfunctionSteps(E)));
    vector<function<Scalar(Scalar, Scalar)>> result;
    auto bases = make_shared<vector<function<ArrayXs(Scalar)>>>(sectorCount);
    for (int i = 0; i < sectorCount; ++i)
        (*bases)[static_cast<size_t>(i)] = sectors[i]->basisCalculator();

    if (!steps->empty()) {
        int cols = (int) steps->at(0).getY(0).cols();
        for (int column = 0; column < cols; ++column) {
            result.push_back([this, steps, E, column, bases](Scalar x, Scalar y) -> Scalar {
                Eigen::Index sectorIndex;
                {
                    Eigen::Index a = 0;
                    Eigen::Index b = this->sectorCount;
                    while (a + 1 < b) {
                        Eigen::Index c = (a + b) / 2;
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