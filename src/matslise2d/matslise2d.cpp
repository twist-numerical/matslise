#include <iostream>
#include <map>
#include <queue>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/find_sector.h"
#include "./matching.h"

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
                               const matslise::Rectangle<2, Scalar> &domain, const Config &_config):
        AbstractMatslise2D<Scalar>(potential, domain), config(_config) {
    dirichletBoundary = Y<Scalar, Eigen::Dynamic>::Dirichlet(config.basisSize);
    auto sectorsBuild = sector_builder::getOrAutomatic(
            _config.ySectorBuilder, _config.tolerance)(this, domain.min, domain.max);
    sectors = std::move(sectorsBuild.sectors);
    matchIndex = sectorsBuild.matchIndex;
    Index sectorCount = sectors.size();
    for (auto &sector : sectors)
        sector->quadratures.reset();

    M.reserve(sectorCount - 1);

    map<pair<int, int>, ArrayXXs> prev;
    map<pair<int, int>, ArrayXXs> next;
    Scalar xWidth = domain.template getMax<0>() - domain.template getMin<0>();
    for (int k = 0; k < sectorCount - 1; ++k) {
        if (k > 0) {
            prev = std::move(next);
            next.clear();
        }

        M.push_back(move(
                gauss_kronrod::adaptive<Scalar, MatrixXs, true>([&, k](const ArrayXs &x) {
                    int depth = static_cast<int>(round(log2(xWidth / (x[x.size() - 1] - x[0]))));
                    int offset = static_cast<int>(round(
                            (x[0] - domain.template getMin<0>()) / (xWidth / (1 << depth))));
                    pair<int, int> key{depth, offset};
                    if (prev.find(key) == prev.end())
                        prev[key] = sectors[k]->template basis<false>(x);
                    const ArrayXXs &prevBasis = prev[key];
                    const ArrayXXs &nextBasis = next[key] = sectors[k + 1]->template basis<false>(x);

                    Array<MatrixXs, Dynamic, 1> result(x.size());
                    for (Index i = 0; i < x.size(); ++i) {
                        result(i) = nextBasis.row(i).matrix().transpose() * prevBasis.row(i).matrix();
                    }
                    return result;
                }, domain.template getMin<0>(), domain.template getMax<0>(), 1e-8, [](const MatrixXs &v) {
                    return v.array().abs().maxCoeff();
                })
        ));
    }
}

template<typename Scalar>
Matslise2D<Scalar>::~Matslise2D() {
    for (auto &sector : sectors)
        delete sector;
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
    Index sectorCount = sectors.size();
    Y<Scalar, Dynamic> yr = sectors[sectorCount - 1]->propagate(
            E, dirichletBoundary, sectors[sectorCount - 1]->max, sectors[sectorCount - 1]->min, use_h);
    conditionY(yr);
    for (int i = sectorCount - 2; i > matchIndex; --i) {
        yr = sectors[i]->propagate(E, (MatrixXs)(M[i].transpose()) * yr, sectors[i]->max, sectors[i]->min, use_h);
        conditionY(yr);
    }

    return errorMatrix<Scalar>(yl, yr);
}


template<typename Scalar>
Index Matslise2D<Scalar>::estimateIndex(Y<Scalar, Eigen::Dynamic> y, const Scalar &E) const {
    Index index = 0;
    Index sectorIndex = 0;
    Index sectorCount = sectors.size();
    for (; sectorIndex < sectorCount - 1; ++sectorIndex) {
        Sector *sector = sectors[sectorIndex];
        Y<Scalar, Dynamic> y1 = sector->propagate(E, y, sector->min, sector->max);
        index += sector->estimateIndex(E, y, y1);
        y = y1;
        conditionY(y);
        y = M[sectorIndex] * y;
    }
    {
        Sector *sector = sectors[sectorIndex];
        index += sector->estimateIndex(E, y, sector->propagate(E, y, sector->min, sector->max));
    }
    return index;
}

template<typename Scalar>
Y<Scalar, Dynamic> Matslise2D<Scalar>::propagate(
        const Scalar &E, const Y<Scalar, Dynamic> &y0, const Scalar &a, const Scalar &b, bool use_h) const {
    if (!domain.contains(1, a) || !domain.contains(1, b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Scalar, Dynamic> y = y0;
    int sectorIndex = find_sector<Matslise2D<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    Index sectorCount = sectors.size();
    while (true) {
        sector = sectors[sectorIndex];
        if (direction == -1 && sectorIndex < sectorCount - 1)
            y = (MatrixXs)(M[sectorIndex].transpose()) * y;
        y = sector->propagate(E, y, a, b, use_h);
        if (sector->contains(b))break;
        conditionY(y);
        if (direction == 1)
            y = M[sectorIndex] * y;
        sectorIndex += direction;
    }
    return y;
}

template<typename Scalar>
vector<pair<Scalar, Scalar>> Matslise2D<Scalar>::matchingErrors(
        const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    pair<MatrixXs, MatrixXs> error_matrix = matchingErrorMatrix(yLeft, E, use_h);
    return eigenvaluesWithDerivatives<Scalar, false>(error_matrix.first, error_matrix.second);
}

template<typename Scalar>
bool newtonSorter(const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) {
    return abs(a.first * b.second) < abs(b.first * a.second);
}

template<typename Scalar>
pair<Scalar, Scalar>
Matslise2D<Scalar>::matchingError(const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    vector<pair<Scalar, Scalar>> errors = matchingErrors(yLeft, E, use_h);
    return *min_element(errors.begin(), errors.end(), &newtonSorter<Scalar>);
}

template<typename Scalar>
Scalar Matslise2D<Scalar>::estimatePotentialMinimum() const {
    auto iterator = this->sectors.begin();
    Scalar minimal = (*iterator++)->matslise->estimatePotentialMinimum();
    for (; iterator != this->sectors.end(); ++iterator)
        minimal = min(minimal, (*iterator)->matslise->estimatePotentialMinimum());
    return minimal;
}

template<typename Scalar>
pair<Scalar, Index> Matslise2D<Scalar>::eigenvalue(const Y<Scalar, Dynamic> &left, const Scalar &_E, bool use_h) const {
    const Scalar tolerance = 1e-9;
    const Scalar minTolerance = 1e-5;
    const int maxIterations = 30;

    Scalar E = _E;
    Scalar error, derror;
    int i = 0;
    do {
        tie(error, derror) = matchingError(left, E, use_h);
        E -= error / derror;
        ++i;
    } while (i < maxIterations && abs(error) > tolerance);

    if (abs(error) > minTolerance)
        return {NAN, 0};

    vector<pair<Scalar, Scalar>> errors = matchingErrors(left, E, use_h);
    sort(errors.begin(), errors.end(), &newtonSorter<Scalar>);
    Index count = 1;
    Index N = config.basisSize;
    for (Index j = 1; j < N; ++j)
        if (abs(errors[j].first / errors[j].second) < minTolerance)
            ++count;
    return {E, count};
}

template<typename Scalar>
Scalar Matslise2D<Scalar>::eigenvalueError(const Y<Scalar, Dynamic> &left, const Scalar &E) const {
    return abs(E - eigenvalue(left, E, false).first);
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> eigenvaluesHelper(
        const Matslise2D<Scalar> &matslise2d, const Y<Scalar, Dynamic> &left,
        const Scalar &Emin, const Scalar &Emax, const Index &Imin, const Index &Imax) {

    queue<tuple<int, Scalar, Scalar, Index, Index>> toCheck;
    toCheck.emplace(0, Emin, Emax, matslise2d.estimateIndex(left, Emin), matslise2d.estimateIndex(left, Emax));
    vector<tuple<Index, Scalar, Index>> found;
    while (!toCheck.empty()) {
        const int &depth = get<0>(toCheck.front());
        const Scalar &a = get<1>(toCheck.front());
        const Scalar &b = get<2>(toCheck.front());
        const Index &ia = get<3>(toCheck.front());
        const Index &ib = get<4>(toCheck.front());

        Scalar mid = (a + b) / 2;
        pair<Scalar, Index> E = matslise2d.eigenvalue(left, mid);
        if (E.second == ib - ia && a - 1e-4 < E.first && E.first < b + 1e-4) {
            if (Imin < ia + E.second || ia < Imax)
                found.emplace_back(ia, E.first, E.second);
        } else {
            Index imid = matslise2d.estimateIndex(left, mid);
            if (imid < ia || imid > ib)
                cerr << "Matslise2D: Error in index estimate" << endl;
            if (imid > ia && imid >= Imin)
                toCheck.emplace(depth + 1, a, mid, ia, imid);
            if (imid < ib && imid <= Imax)
                toCheck.emplace(depth + 1, mid, b, imid, ib);
        }
        toCheck.pop();
    }
    sort(found.begin(), found.end());
    return found;
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> Matslise2D<Scalar>::eigenvalues(
        const Y<Scalar, Dynamic> &left, const Scalar &Emin, const Scalar &Emax) const {
    return eigenvaluesHelper(*this, left, Emin, Emax, 0, numeric_limits<Index>::max());
}

template<typename Scalar>
vector<tuple<Index, Scalar, Index>> Matslise2D<Scalar>::eigenvaluesByIndex(
        const Y<Scalar, Dynamic> &left, Index Imin, Index Imax) const {
    if (Imin < 0)
        Imin = 0;
    if (Imax <= Imin)
        return vector<tuple<Index, Scalar, Index>>();
    Scalar step = 1;
    Scalar Emin = estimatePotentialMinimum();
    Index i = estimateIndex(left, Emin);
    while (i > Imin) {
        Emin -= step;
        step *= 2;
        i = estimateIndex(left, Emin);
    }
    step = 1;
    Scalar Emax = Emin;
    do {
        Emax += step;
        step *= 2;
        i = estimateIndex(left, Emax);
        if (i < Imin)
            Emin = Emax;
    } while (i < Imax);

    return eigenvaluesHelper(*this, left, Emin, Emax, Imin, Imax);
}

#include "../util/instantiate.h"
