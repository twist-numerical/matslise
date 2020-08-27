#include <complex>
#include <queue>
#include "../matslise.h"
#include "../util/constants.h"
#include "../util/calculateEta.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar, typename Derived>
Matrix<complex<Scalar>, Dynamic, Dynamic> theta(
        const MatrixBase<Derived> &U,
        const MatrixBase<Derived> &V) {
    MATSLISE_SCOPED_TIMER("ND ios: theta");
    return (V - U * complex<Scalar>(0, 1))
            .transpose()
            .partialPivLu()
            .solve((V + U * complex<Scalar>(0, 1)).transpose())
            .transpose();
}

template<typename Scalar>
inline Array<Scalar, Dynamic, 1> angle(const Matrix<complex<Scalar>, Dynamic, Dynamic> &m) {
    MATSLISE_SCOPED_TIMER("ND ios: angle");
    const Scalar PI2 = constants<Scalar>::PI * 2;
    return m.eigenvalues().array().arg().unaryExpr(
            [&](const Scalar &a) -> Scalar { return a < -1e-16 ? a + PI2 : a; });
}

template<typename Scalar>
Index estimateIndexOfSector(const typename Matscs<Scalar>::Sector &sector,
                            const Scalar &E, const Y<Scalar, Eigen::Dynamic> &y0, const Y<Scalar, Eigen::Dynamic> &y1) {
    MATSLISE_SCOPED_TIMER("ND index of sector");
    Index n = sector.n;
    using ArrayXs = typename Matslise2D<Scalar>::ArrayXs;
    using MatrixXs = Matrix<Scalar, Dynamic, Dynamic>;
    using MatrixXcs = Matrix<complex<Scalar>, Dynamic, Dynamic>;

    Index zeros = 0;
    queue<tuple<int, Scalar, Scalar, MatrixXs, MatrixXs, ArrayXs>> todo;
    todo.emplace(0, sector.min, sector.max, y0.getY(0), y0.getY(1),
                 angle<Scalar>(theta<Scalar>(y1.getY(0), y1.getY(1))));

    while (!todo.empty()) {
        const int &depth = get<0>(todo.front());
        const Scalar &a = get<1>(todo.front());
        const Scalar &b = get<2>(todo.front());
        Scalar h = b - a;
        const MatrixXs &U0 = get<3>(todo.front());
        const MatrixXs &V0 = get<4>(todo.front());
        const ArrayXs &betas = get<5>(todo.front());

        ArrayXs Z = h * h * (sector.vs[0].diagonal().array() - ArrayXs::Constant(n, E));
        Array<Scalar, 2, Dynamic> eta = calculateEta<Scalar, 2>(Z);

        MatrixXcs thetaZ0 = theta<Scalar>(U0, V0);
        const complex<Scalar> i_delta(0, h);
        ArrayXs alphas = angle<Scalar>(
                ((eta.row(0) + i_delta * eta.row(1)) / (eta.row(0) - i_delta * eta.row(1))).matrix().asDiagonal() *
                thetaZ0
        );
        if (depth < 3 && ((betas < 1e-1).any() || ((betas - alphas).abs() > 5.5).any())) {
            // Zeroth order propagation probably inaccurate: refine steps
            Scalar mid = (a + b) / 2;
            Y<Scalar, Dynamic> yMid = sector.direction == matslise::forward
                                      ? sector.propagateColumn(E, y0, sector.min, mid)
                                      : sector.propagateColumn(E, y1, sector.max, mid);
            todo.emplace(depth + 1, a, mid, U0, V0, angle<Scalar>(theta<Scalar>(yMid.getY(0), yMid.getY(1))));
            todo.emplace(depth + 1, mid, b, yMid.getY(0), yMid.getY(1), betas);
        } else {
            Scalar argdet = 0;
            for (int i = 0; i < n; ++i) {
                if (Z[i] < 0) {
                    Scalar sZ = sqrt(-Z[i]);
                    argdet += (sZ + atan2(
                            (h - sZ) * eta(1, i) * eta(0, i),
                            1 + (h * sZ + Z[i]) * eta(1, i) * eta(1, i)));
                } else {
                    argdet += atan2(h * eta(1, i), eta(0, i));
                }
            }
            argdet *= 2;

            zeros += (Eigen::Index) round(
                    (angle<Scalar>(thetaZ0).sum() + argdet - alphas.sum()) / (2 * constants<Scalar>::PI));
        }
        todo.pop();
    }
    return zeros;
}

template<typename Scalar>
pair<Y<Scalar, Eigen::Dynamic>, Index> MatsliseNDSector<Scalar>::propagateWithIndex(
        const Scalar &E, Y<Scalar, Eigen::Dynamic> y0) const {
    Index index = 0;
    for (auto &sector : matscs) {
        Y<Scalar, Eigen::Dynamic> y1 = sector.propagateColumn(E, y0, min, max, true);
        index += estimateIndexOfSector(sector, E, y0, y1);
        y0 = y1;
    }
    return {y0, index};
}

template<typename Scalar>
template<int r>
Y<Scalar, Eigen::Dynamic, r>
MatsliseNDSector<Scalar>::propagate(
        const Scalar &E, Y<Scalar, Eigen::Dynamic, r> y, const Scalar &a, const Scalar &b, bool use_h) const {
    MATSLISE_SCOPED_TIMER("ND sector propagate");
    if (a < b)
        for (auto sector = matscs.begin(); sector != matscs.end(); ++sector)
            y = sector->propagateColumn(E, y, a, b, use_h);
    else if (b < a)
        for (auto sector = matscs.rbegin(); sector != matscs.rend(); ++sector)
            y = sector->propagateColumn(E, y, a, b, use_h);
    return y;
}

template<typename Scalar>
Scalar MatsliseNDSector<Scalar>::error() const {
    Scalar error = 0;
    for (auto &sector : matscs) {
        error += sector.error();
    }
    return error;
}


#define INSTANTIATE_PROPAGATE(Scalar, r) \
template Y<Scalar, Dynamic, r> \
MatsliseNDSector<Scalar>::propagate<r>( \
        const Scalar &E, Y<Scalar, Eigen::Dynamic, r>, const Scalar &, const Scalar &, bool) const;

#define INSTANTIATE_MORE(Scalar)\
INSTANTIATE_PROPAGATE(Scalar, 1)\
INSTANTIATE_PROPAGATE(Scalar, -1)

#include "../util/instantiate.h"
