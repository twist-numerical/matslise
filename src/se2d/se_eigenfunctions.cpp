#include "../se2d.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;

template<int n>
Y<MatrixXd> *SEBase<n>::computeEigenfunctionSteps(double E) const {
    int match = sectorCount / 2;

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

    MatrixXd big = MatrixXd::Zero(2 * N, 2 * N);
    big << matchLeft.y[0], -matchRight.y[0], matchLeft.y[1], -matchRight.y[1];

    FullPivLU<MatrixXd> lu(big);
    MatrixXd kernel = lu.kernel();

    if (kernel.isZero(0)) {
        delete[] steps;
        return nullptr;
    }

    MatrixXd left = kernel.topRows(N);
    MatrixXd right = kernel.bottomRows(N);
    for (int i = 0; i <= match; ++i)
        steps[i] *= left;
    for (int i = sectorCount; i > match; --i)
        steps[i] *= right;

    return steps;
};


std::vector<typename dim<2>::array>
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


    vector<ArrayXXd> result;
    Y<MatrixXd> *steps = computeEigenfunctionSteps(E);
    if (steps != nullptr) {
        int cols = (int) steps[0].y[0].cols();
        for (int i = 0; i < cols; ++i)
            result.push_back(ArrayXXd::Zero(nx, ny));

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
                    result[i].col(nextY) = prod.col(i);
                ++nextY;
            }
        }

        delete[] steps;
    }

    return result;
}