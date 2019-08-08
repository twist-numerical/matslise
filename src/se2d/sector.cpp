#include <iostream>
#include "../matslise.h"
#include "../util/lobatto.h"

using namespace matslise;
using namespace std;
using namespace Eigen;
using namespace matslise::SEnD_util;
using namespace matslise::sectorbuilder;

template<typename Scalar>
SE2D<Scalar>::Sector::Sector(SE2D<Scalar> *se2d, const Scalar &ymin, const Scalar &ymax, bool backward)
        : se2d(se2d), min(ymin), max(ymax) {
    const Y<Scalar> y0 = Y<Scalar>({0, 1}, {0, 0});

    const Scalar ybar = (ymax + ymin) / 2;
    function<Scalar(Scalar)> vbar_fun = [se2d, ybar](Scalar x) -> Scalar { return se2d->V(x, ybar); };
    vbar = se2d->grid.unaryExpr(vbar_fun);
    if (se2d->options.nestedOptions._symmetric)
        matslise = new HalfRange<Scalar>(vbar_fun, se2d->domain.sub.max, se2d->options.nestedOptions._builder);
    else
        matslise = new Matslise<Scalar>(vbar_fun, se2d->domain.sub.min, se2d->domain.sub.max,
                                        se2d->options.nestedOptions._builder);

    vector<pair<int, Scalar>> index_eigv = matslise->computeEigenvaluesByIndex(0, se2d->N, y0);
    if (static_cast<int>(index_eigv.size()) != se2d->N) {
        throw std::runtime_error("SE2D: not enough basis-functions found on a sector");
    }
    eigenvalues = new Scalar[se2d->N];
    eigenfunctions = new ArrayXs[se2d->N];
    Scalar E;
    int index;
    for (int i = 0; i < se2d->N; ++i) {
        tie(index, E) = index_eigv[static_cast<unsigned long>(i)];
        eigenvalues[i] = E;
        // TODO: check i == index
        Array<Y<Scalar>, Dynamic, 1> func = matslise->computeEigenfunction(E, y0, se2d->grid, index);
        eigenfunctions[i] = ArrayXs(func.size());
        for (int j = 0; j < func.size(); ++j)
            eigenfunctions[i][j] = func[j].y[0];
    }

    matscs = new Matscs<Scalar>(
            [this](Scalar y) -> MatrixXs { return this->calculateDeltaV(y); },
            se2d->N, ymin, ymax,
            std::shared_ptr<SectorBuilder<Matscs<Scalar>>>(
                    new Custom<Matscs<Scalar>>([this, backward](Matscs<Scalar> *p, Scalar xmin, Scalar xmax) {
                        int n = this->se2d->options._stepsPerSector;
                        Scalar h = (xmax - xmin) / n;
                        p->sectorCount = n;
                        p->sectors = new typename Matscs<Scalar>::Sector *[n];
                        Scalar left = xmin;
                        for (int i = 0; i < n; ++i) {
                            Scalar right = xmax - (n - i - 1) * h;
                            p->sectors[i] = new typename Matscs<Scalar>::Sector(p, left, right, backward);
                            left = right;
                        }
                        p->match = p->sectors[n - 1]->min;
                    })));
}

template<typename Scalar>
SE2D<Scalar>::Sector::~Sector() {
    delete matslise;
    delete matscs;
    delete[] eigenvalues;
    delete[] eigenfunctions;
}

template<typename Scalar>
template<int r>
Y<Scalar, Eigen::Dynamic, r>
SE2D<Scalar>::Sector::propagate(
        const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b, bool use_h) const {
    return matscs->propagateColumn(
            E, y0, a < min ? min : a > max ? max : a, b < min ? min : b > max ? max : b, use_h);
}

template<typename Scalar>
Scalar SE2D<Scalar>::Sector::calculateError() const {
    Scalar error = 0;
    for (int i = 0; i < matscs->sectorCount; ++i)
        error += matscs->sectors[i]->calculateError();
    return error;
}

template<typename Scalar>
typename SE2D<Scalar>::MatrixXs SE2D<Scalar>::Sector::calculateDeltaV(const Scalar &y) const {
    MatrixXs dV(se2d->N, se2d->N);

    ArrayXs vDiff = se2d->grid.unaryExpr([this, y](const Scalar &x) -> Scalar { return this->se2d->V(x, y); }) - vbar;

    for (int i = 0; i < se2d->N; ++i) {
        for (int j = 0; j <= i; ++j) {
            dV(i, j) = lobatto::quadrature<Scalar>(se2d->grid, eigenfunctions[i] * vDiff * eigenfunctions[j]);
            if (j < i) dV(j, i) = dV(i, j);
        }
        dV(i, i) += eigenvalues[i];
    }

    return dV;
}

template<typename Scalar>
typename SE2D<Scalar>::ArrayXs
SE2D<Scalar>::Sector::computeEigenfunction(int index, const typename SE2D<Scalar>::ArrayXs &x) const {
    const Y<Scalar> y0 = Y<Scalar>({0, 1}, {0, 0});
    Eigen::Index size = x.size();

    Array<matslise::Y<Scalar>, Dynamic, 1> raw = matslise->computeEigenfunction(eigenvalues[index], y0, x, index);
    ArrayXs result(size);
    for (int i = 0; i < size; ++i)
        result(i) = raw(i).y[0];
    return result;
}

template<typename Scalar>
function<typename SE2D<Scalar>::ArrayXs(Scalar)> SE2D<Scalar>::Sector::basisCalculator() const {
    const Y<Scalar> y0 = Y<Scalar>::Dirichlet(1);
    vector<function<Y<Scalar>(Scalar)>> basis(static_cast<size_t>(se2d->N));
    for (int index = 0; index < se2d->N; ++index) {
        basis[static_cast<size_t>(index)] = matslise->eigenfunctionCalculator(eigenvalues[index], y0, index);
    }
    return [basis](const Scalar &x) -> ArrayXs {
        ArrayXs result(basis.size());
        for (int i = 0; i < static_cast<int>(basis.size()); ++i)
            result[i] = basis[static_cast<size_t>(i)](x).y(0, 0);
        return result;
    };
}

#define INSTANTIATE_PROPAGATE(Scalar, r)\
template Y<Scalar, Dynamic, r>\
SE2D<Scalar>::Sector::propagate<r>(\
        const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b, bool use_h) const;

#define INSTANTIATE_MORE(Scalar)\
INSTANTIATE_PROPAGATE(Scalar, 1)\
INSTANTIATE_PROPAGATE(Scalar, -1)

#include "../util/instantiate.h"