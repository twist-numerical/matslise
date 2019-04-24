#include "../se2d.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;

#define cec_cce(y) ((y).getdY(0).transpose()*(y).getY(1) - (y).getY(0).transpose()*(y).getdY(1))

template<int n>
Y<Dynamic> *SEnD<n>::computeEigenfunctionSteps(double E) const {
    Y<Dynamic> *steps = new Y<Dynamic>[sectorCount + 1];

    steps[0] = Y<Dynamic>::Dirichlet(N);
    steps[sectorCount] = Y<Dynamic>::Dirichlet(N);

    MatrixXd normRight = MatrixXd::Zero(N, N);
    int matchIndex = 0;
    for (int i = sectorCount - 1; sectors[i]->min > match; --i) {
        Y<Dynamic> next = i < sectorCount - 1 ? (MatrixXd) (M[i].transpose()) * steps[i + 1] : steps[i + 1];
        steps[i] = sectors[i]->propagate(E, next, false);
        normRight += cec_cce(next) - cec_cce(steps[i]);
        matchIndex = i;
    }
    Y<Dynamic> matchRight = steps[matchIndex];

    MatrixXd normLeft = MatrixXd::Zero(N, N);
    for (int i = 0; i < matchIndex; ++i) {
        Y<Dynamic> next = sectors[i]->propagate(E, steps[i], true);
        steps[i + 1] = M[i] * next;
        normLeft += cec_cce(next) - cec_cce(steps[i]);
    }
    Y<Dynamic> matchLeft = steps[matchIndex];

    ColPivHouseholderQR<MatrixXd> left_solver(matchLeft.getY(0).transpose());
    ColPivHouseholderQR<MatrixXd> right_solver(matchRight.getY(0).transpose());
    MatrixXd Ul = left_solver.solve(matchLeft.getY(1).transpose()).transpose();
    MatrixXd Ur = right_solver.solve(matchRight.getY(1).transpose()).transpose();

    FullPivLU<MatrixXd> lu(Ul - Ur);
    lu.setThreshold(1e-4);
    if (lu.dimensionOfKernel() == 0) {
        delete[] steps;
        return nullptr;
    }
    MatrixXd kernel = lu.kernel();

    Y<Dynamic> *elements = new Y<Dynamic>[sectorCount + 1];
    MatrixXd left = matchLeft.getY(0).colPivHouseholderQr().solve(kernel);
    MatrixXd right = matchRight.getY(0).colPivHouseholderQr().solve(kernel);
    VectorXd scaling = (left.transpose() * normLeft * left
                        + right.transpose() * normRight * right).diagonal();
    scaling = scaling.unaryExpr([](double s) { return s < 0 ? 1 : 1. / sqrt(s); });
    left *= scaling.asDiagonal();
    right *= scaling.asDiagonal();
    for (int i = 0; i <= matchIndex; ++i)
        elements[i] = steps[i] * left;
    for (int i = sectorCount; i > matchIndex; --i)
        elements[i] = steps[i] * right;
    delete[] steps;

    return elements;
};

template<>
std::vector<typename dim<2>::array> *
SEnD<2>::computeEigenfunction(double E, const Eigen::ArrayXd (&xs)[2]) const {
    const Eigen::ArrayXd &x = xs[0];
    const Eigen::ArrayXd &y = xs[1];
    long nx = x.size();
    for (int i = 1; i < nx; ++i)
        if (x[i - 1] > x[i])
            throw runtime_error("SE2D::computeEigenfunction(): x has to be sorted");

    long ny = y.size();
    for (int i = 1; i < ny; ++i)
        if (y[i - 1] > y[i])
            throw runtime_error("SE2D::computeEigenfunction(): y has to be sorted");

    if (y[0] < sectors[0]->min || y[ny - 1] > sectors[sectorCount - 1]->max)
        throw runtime_error("SE2D::computeEigenfunction(): y is out of range");


    auto result = new vector<ArrayXXd>;
    Y<Dynamic> *steps = computeEigenfunctionSteps(E);
    if (steps != nullptr) {
        int cols = (int) steps[0].getY(0).cols();
        for (int i = 0; i < cols; ++i)
            result->push_back(ArrayXXd::Zero(nx, ny));

        int nextY = 0;
        int sector = 0;
        while (nextY < ny) {
            while (sector < sectorCount && y[nextY] > sectors[sector]->max) {
                ++sector;
                if (sector >= sectorCount)
                    throw runtime_error("SE2D::computeEigenfunction(): y is out of range");
            }

            MatrixXd B(nx, N);
            for (int j = 0; j < N; ++j)
                B.col(j) = sectors[sector]->computeEigenfunction(j, x);

            while (nextY < ny && y[nextY] <= sectors[sector]->max) {
                MatrixXd prod = B * sectors[sector]->propagate(E, steps[sector], y[nextY], true).getY(0);
                for (int i = 0; i < cols; ++i)
                    result->at(i).col(nextY) = prod.col(i);
                ++nextY;
            }
        }

        delete[] steps;
    }

    return result;
}