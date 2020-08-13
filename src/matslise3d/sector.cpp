#include "../matslise.h"
#include "../util/legendre.h"
#include "../util/quadrature.h"

using namespace Eigen;
using namespace matslise;
using namespace std;
using namespace quadrature;

template<typename Scalar>
Array<Scalar, Dynamic, Dynamic>
eval2d(const function<Scalar(const Scalar &, const Scalar &)> &f,
       const Array<Scalar, Dynamic, 1> &x, const Array<Scalar, Dynamic, 1> &y) {
    Array<Scalar, Dynamic, Dynamic> result(x.size(), y.size());
    for (Index i = 0; i < x.size(); ++i)
        for (Index j = 0; j < x.size(); ++j)
            result(i, j) = f(x[i], y[j]);
    return result;
}

template<typename Scalar>
vector<typename Matscs<Scalar>::Sector> initializeMatscs(const typename Matslise3D<Scalar>::Sector &sector) {
    const Index &N = sector.matslise3d->config.xyBasisSize;
    auto &matslise3d = sector.matslise3d;
    typedef typename Matslise3D<Scalar>::MatrixXs MatrixXs;
    typedef typename Matslise3D<Scalar>::ArrayXXs ArrayXXs;

    vector<typename Matscs<Scalar>::Sector> matscs;
    Index steps = matslise3d->config.zStepsPerSector;
    Scalar h = (sector.max - sector.min) / steps;
    matscs.reserve(steps);
    for (Index i = 0; i < steps; ++i) {
        Scalar min = sector.min + i * h;
        Scalar max = sector.max - (steps - i - 1) * h;
        // cout << min << ", " << max << endl; To many calls
        matscs.emplace_back(
                legendre::getCoefficients<MATSCS_N, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Scalar>(
                        [&](Scalar z) -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> {
                            MatrixXs dV(N, N);

                            ArrayXXs vDiff = eval2d<Scalar>([&](const Scalar &x, const Scalar &y) {
                                return matslise3d->potential(x, y, z);
                            }, matslise3d->grid_x, matslise3d->grid_y) - sector.vbar;

                            for (int i = 0; i < N; ++i) {
                                for (int j = 0; j <= i; ++j) {
                                    dV(i, j) = lobatto::quadrature<Scalar>(
                                            matslise3d->grid_x, matslise3d->grid_y,
                                            sector.eigenfunctions_grid[i] * vDiff * sector.eigenfunctions_grid[j]);
                                    if (j < i) dV(j, i) = dV(i, j);
                                }
                                dV(i, i) += sector.eigenvalues[i];
                            }

                            // cout << "Z: " << z << " dV:" << endl;
                            // cout << dV << "\n" << endl;
                            return dV;
                        }, min, max),
                min, max, sector.direction);
    }
    return matscs;
}

template<typename Scalar>
Matslise3D<Scalar>::Sector::Sector(
        const Matslise3D<Scalar> *matslise3d, const Scalar &zmin, const Scalar &zmax, Direction direction)
        : matslise3d(matslise3d), min(zmin), max(zmax), direction(direction) {

    zbar = (zmax + zmin) / 2;
    function<Scalar(const Scalar &, const Scalar &)> vbar_fun = [matslise3d, this](
            const Scalar &x, const Scalar &y) -> Scalar {
        return matslise3d->potential(x, y, zbar);
    };

    matslise2d = std::make_shared<Matslise2D<Scalar>>(
            vbar_fun, matslise3d->domain.template slice<0, 1>(),
            [&] {
                typename Matslise2D<Scalar>::Config config2d;
                config2d.tolerance = matslise3d->config.tolerance;
                config2d.basisSize = matslise3d->config.xBasisSize;
                config2d.stepsPerSector = matslise3d->config.yStepsPerSector;
                return config2d;
            }());

    const Index &N = matslise3d->config.xyBasisSize;
    cout << "Seeking: " << zbar << endl;
    vector<tuple<Index, Scalar, Index>> singleEigenvalues = matslise2d->eigenvaluesByIndex(0, N);
    cout << "found: " << zbar << endl;
    Index eigenvalueCount = 0;
    for (auto &iEm : singleEigenvalues) {
        eigenvalueCount += get<2>(iEm);
        cout << " (" << get<0>(iEm) << ", " << get<1>(iEm) << ", " << get<2>(iEm) << ")";
    }
    cout << endl;
    if (eigenvalueCount < N) {
        throw std::runtime_error("Matlise3D: not enough basis-functions found on a sector");
    }
    eigenvalues.reserve(N);
    eigenfunctions.reserve(N);
    eigenfunctions_grid.reserve(N);

    Index i = 0;
    for (auto &E : singleEigenvalues) {
        for (auto &f :  matslise2d->eigenfunction(get<1>(E))) {
            eigenvalues.push_back(get<1>(E));
            eigenfunctions.push_back(f);
            eigenfunctions_grid.push_back(move(f(matslise3d->grid_x, matslise3d->grid_y)));

            if (++i == N)
                break;
        }
        if (i == N)
            break;
    }
    vbar = eval2d<Scalar>([&](const Scalar &x, const Scalar &y) {
        return matslise3d->potential(x, y, zbar);
    }, matslise3d->grid_x, matslise3d->grid_y);

    matscs = initializeMatscs<Scalar>(*this);
}

template<typename Scalar>
void Matslise3D<Scalar>::Sector::setDirection(Direction newDirection) {
    direction = newDirection;
    for (auto &sector : matscs)
        sector.setDirection(newDirection);
}

template<typename Scalar>
template<int r>
Y<Scalar, Eigen::Dynamic, r>
Matslise3D<Scalar>::Sector::propagate(
        const Scalar &E, Y<Scalar, Eigen::Dynamic, r> y, const Scalar &a, const Scalar &b, bool use_h) const {
    if (a < b)
        for (auto sector = matscs.begin(); sector != matscs.end(); ++sector)
            y = sector->propagateColumn(E, y, a, b, use_h);
    else if (b < a)
        for (auto sector = matscs.rbegin(); sector != matscs.rend(); ++sector)
            y = sector->propagateColumn(E, y, a, b, use_h);
    return y;
}

template<typename Scalar>
Scalar Matslise3D<Scalar>::Sector::error() const {
    Scalar error = 0;
    for (auto &sector : matscs)
        error += sector.error();
    return error;
}


template<typename Scalar>
typename Matslise3D<Scalar>::Sector *Matslise3D<Scalar>::Sector::refine(
        const Matslise3D<Scalar> *problem, const Scalar &_min, const Scalar &_max, Direction _direction) const {
    Scalar h = _max - _min;
    if (zbar < _min + h / 3 || zbar > _max - h / 3) {
        return new Sector(problem, _min, _max, _direction);
    }

    // std::cout << "refining sector 2d" << std::endl;
    auto sector = new Sector(problem);
    sector->direction = _direction;
    sector->min = _min;
    sector->max = _max;
    sector->zbar = zbar;
    sector->matslise2d = matslise2d;
    sector->eigenfunctions = eigenfunctions;
    sector->eigenvalues = eigenvalues;
    sector->eigenfunctions_grid = eigenfunctions_grid;
    sector->matscs = initializeMatscs<Scalar>(*sector);
    return sector;
}

#define INSTANTIATE_PROPAGATE(Scalar, r) \
template Y<Scalar, Dynamic, r> \
Matslise3D<Scalar>::Sector::propagate<r>( \
        const Scalar &, Y<Scalar, Eigen::Dynamic, r>, const Scalar &, const Scalar &, bool) const;

#define INSTANTIATE_MORE(Scalar) \
INSTANTIATE_PROPAGATE(Scalar, 1) \
INSTANTIATE_PROPAGATE(Scalar, -1)

#include "../util/instantiate.h"