#include "../se2d.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;

template<int n>
Y<Dynamic> *SEBase<n>::computeEigenfunctionSteps(double E) const {
    Y<Dynamic> *steps = new Y<Dynamic>[sectorCount + 1];

    steps[0] = Y<Dynamic>::Dirichlet(N);
    steps[sectorCount] = Y<Dynamic>::Dirichlet(N);

    for (int i = sectorCount - 1; i >= match; --i)
        steps[i] = sectors[i]->propagate(
                E, i < sectorCount - 1 ? (MatrixXd) (M[i].transpose()) * steps[i + 1] : steps[i + 1], false);
    Y<Dynamic> matchRight = steps[match];

    for (int i = 0; i < match; ++i)
        steps[i + 1] = M[i] * sectors[i]->propagate(E, steps[i], true);
    Y<Dynamic> matchLeft = steps[match];

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
    for (int i = 0; i <= match; ++i)
        elements[i] = steps[i] * left;
    for (int i = sectorCount; i > match; --i)
        elements[i] = steps[i] * right;
    delete[] steps;

    return elements;
};


std::vector<typename dim<2>::array> *
SEnD<2>::computeEigenfunction(double E, const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const {
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