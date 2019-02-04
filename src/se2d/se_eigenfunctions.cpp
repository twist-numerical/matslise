#include "../se2d.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;

template<int n>
Y<MatrixXd> *SEBase<n>::computeEigenfunctionSteps(double E) const {
    int match = sectorCount / 2 + 1;

    Y<MatrixXd> *steps = new Y<MatrixXd>[sectorCount + 1];

    steps[0] = Y<MatrixXd>({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)},
                           {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});
    steps[sectorCount] = Y<MatrixXd>({MatrixXd::Zero(N, N), MatrixXd::Identity(N, N)},
                                     {MatrixXd::Zero(N, N), MatrixXd::Zero(N, N)});

    for (int i = sectorCount - 1; i >= match; --i)
        steps[i] = sectors[i]->propagate(
                E, i < sectorCount - 1 ? M[i].transpose() * steps[i + 1] : steps[i + 1], false);
    Y<MatrixXd> matchRight = steps[match];

    for (int i = 0; i < match; ++i)
        steps[i + 1] = M[i] * sectors[i]->propagate(E, steps[i], true);
    Y<MatrixXd> matchLeft = steps[match];

    ColPivHouseholderQR<MatrixXd> left_solver(matchLeft.y[0].transpose());
    ColPivHouseholderQR<MatrixXd> right_solver(matchRight.y[0].transpose());
    MatrixXd Ul = left_solver.solve(matchLeft.y[1].transpose()).transpose();
    MatrixXd Ur = right_solver.solve(matchRight.y[1].transpose()).transpose();

    FullPivLU<MatrixXd> lu(Ul - Ur);
    lu.setThreshold(1e-4);
    if (lu.dimensionOfKernel() == 0) {
        delete[] steps;
        return nullptr;
    }
    MatrixXd kernel = lu.kernel();

    MatrixXd left = matchLeft.y[0].colPivHouseholderQr().solve(kernel);
    MatrixXd right = matchRight.y[0].colPivHouseholderQr().solve(kernel);
    for (int i = 0; i <= match; ++i)
        steps[i] *= left;
    for (int i = sectorCount; i > match; --i)
        steps[i] *= right;

    return steps;
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
    Y<MatrixXd> *steps = computeEigenfunctionSteps(E);
    if (steps != nullptr) {
        int cols = (int) steps[0].y[0].cols();
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
                MatrixXd prod = B * sectors[sector]->propagate(E, steps[sector], y[nextY], true).y[0];
                for (int i = 0; i < cols; ++i)
                    result->at(i).col(nextY) = prod.col(i);
                ++nextY;
            }
        }

        delete[] steps;
    }

    return result;
}