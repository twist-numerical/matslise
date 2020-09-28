#include <complex>
#include <queue>
#include "../matslise.h"
#include "../util/constants.h"
#include "../util/calculateEta.h"
#include "../util/find_sector.h"
#include "../util/matching.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic> conditionY(Y<Scalar, Dynamic> &y) {
    MATSLISE_SCOPED_TIMER("ND conditionY");
    Matrix<Scalar, Dynamic, Dynamic> U = y.getY(0).partialPivLu().matrixLU();
    U.template triangularView<StrictlyLower>().setZero();
    U.template triangularView<Upper>().
            template solveInPlace<OnTheRight>(y.data);
    return U;
}

template<typename Scalar, typename Sector>
Y<Scalar, Dynamic> MatsliseND<Scalar, Sector>::propagate(
        const Scalar &E, const Y<Scalar, Dynamic> &y0, const Scalar &a, const Scalar &b, bool use_h) const {
    MATSLISE_SCOPED_TIMER("ND propagate");
    Y<Scalar, Dynamic> y = y0;
    int sectorIndex = find_sector<MatsliseND<Scalar, Sector>>(this, a);
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


template<typename Scalar, typename Sector>
pair<typename MatsliseND<Scalar, Sector>::MatrixXs, typename MatsliseND<Scalar, Sector>::MatrixXs>
MatsliseND<Scalar, Sector>::matchingErrorMatrix(
        const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    MATSLISE_SCOPED_TIMER("ND matchingErrorMatrix");
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


template<typename Scalar, typename Sector>
vector<pair<Scalar, Scalar>> MatsliseND<Scalar, Sector>::matchingErrors(
        const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    MATSLISE_SCOPED_TIMER("ND matchingErrors");
    pair<MatrixXs, MatrixXs> error_matrix = matchingErrorMatrix(yLeft, E, use_h);
    return eigenvaluesWithDerivatives<Scalar, false>(error_matrix.first, error_matrix.second);
}

template<typename Scalar>
bool newtonSorter(const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) {
    return abs(a.first * b.second) < abs(b.first * a.second);
}

template<typename Scalar, typename Sector>
pair<Scalar, Scalar>
MatsliseND<Scalar, Sector>::matchingError(const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    vector<pair<Scalar, Scalar>> errors = matchingErrors(yLeft, E, use_h);
    return *min_element(errors.begin(), errors.end(), &newtonSorter<Scalar>);
}

template<typename Scalar, typename Sector>
pair<Scalar, Index> MatsliseND<Scalar, Sector>::eigenvalue(
        const Y<Scalar, Dynamic> &left, const Scalar &_E, bool use_h) const {
    MATSLISE_SCOPED_TIMER("ND eigenvalue");
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
    for (Index j = 1; j < basisSize; ++j)
        if (abs(errors[j].first / errors[j].second) < minTolerance)
            ++count;
    return {E, count};
}

template<typename Scalar, typename Sector>
Scalar MatsliseND<Scalar, Sector>::eigenvalueError(const Y<Scalar, Dynamic> &left, const Scalar &E) const {
    return abs(E - eigenvalue(left, E, false).first);
}

template<typename Scalar, typename Sector>
vector<tuple<Index, Scalar, Index>> eigenvaluesHelper(
        const MatsliseND<Scalar, Sector> &matslise2d, const Y<Scalar, Dynamic> &left,
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
        } else if (depth > 30) {
            cerr << "MatsliseND: max search depth reached" << endl;
        } else {
            Index imid = matslise2d.estimateIndex(left, mid);
            if (imid < ia || imid > ib)
                cerr << "MatsliseND: Error in index estimate" << endl;
            if (imid > ia && imid > Imin)
                toCheck.emplace(depth + 1, a, mid, ia, imid);
            if (imid < ib && imid <= Imax)
                toCheck.emplace(depth + 1, mid, b, imid, ib);
        }
        toCheck.pop();
    }
    sort(found.begin(), found.end());
    return found;
}

template<typename Scalar, typename Sector>
vector<tuple<Index, Scalar, Index>> MatsliseND<Scalar, Sector>::eigenvalues(
        const Y<Scalar, Dynamic> &left, const Scalar &Emin, const Scalar &Emax) const {
    MATSLISE_SCOPED_TIMER("ND eigenvalues");
    return eigenvaluesHelper(*this, left, Emin, Emax, 0, numeric_limits<Index>::max());
}

template<typename Scalar, typename Sector>
vector<tuple<Index, Scalar, Index>> MatsliseND<Scalar, Sector>::eigenvaluesByIndex(
        const Y<Scalar, Dynamic> &left, Index Imin, Index Imax) const {
    MATSLISE_SCOPED_TIMER("ND eigenvaluesByIndex");
    if (Imin < 0)
        Imin = 0;
    if (Imax <= Imin)
        return vector<tuple<Index, Scalar, Index>>();
    Scalar step = 1;
    Scalar Emin = 0;
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

template<typename Scalar, typename Sector>
Index MatsliseND<Scalar, Sector>::estimateIndex(Y<Scalar, Eigen::Dynamic> y, const Scalar &E) const {
    MATSLISE_SCOPED_TIMER("ND estimateIndex");
    Index index = 0;
    Index sectorIndex = 0;
    Index tmpIndex;
    Index sectorCount = sectors.size();
    for (; sectorIndex < sectorCount - 1; ++sectorIndex) {
        Sector *sector = sectors[sectorIndex];
        tie(y, tmpIndex) = sector->propagateWithIndex(E, y);
        index += tmpIndex;
        conditionY(y);
        y = M[sectorIndex] * y;
    }
    index += sectors[sectorIndex]->propagateWithIndex(E, y).second;
    return index;
}

template<typename Scalar, bool add>
void cec_cce(Matrix<Scalar, Dynamic, Dynamic> &addTo, const Y<Scalar, Dynamic, Dynamic> &y) {
    if constexpr(add)
        addTo += (y).getdY(0).transpose() * (y).getY(1) - (y).getdY(1).transpose() * (y).getY(0);
    else
        addTo -= (y).getdY(0).transpose() * (y).getY(1) - (y).getdY(1).transpose() * (y).getY(0);
}

template<typename Scalar, typename Sector>
vector<Y<Scalar, Dynamic>>
MatsliseND<Scalar, Sector>::eigenfunctionSteps(const Y<Scalar, Dynamic> &yLeft, const Scalar &E) const {
    MATSLISE_SCOPED_TIMER("ND eigenfunctionSteps");
    Index sectorCount = sectors.size();
    auto *steps = new Y<Scalar, Dynamic>[sectorCount + 1];
    auto *endSteps = new Y<Scalar, Dynamic>[sectorCount];

    steps[0] = yLeft;
    steps[sectorCount] = dirichletBoundary;
    auto *U = new MatrixXs[sectorCount];

    for (int i = sectorCount - 1; i > matchIndex; --i) {
        endSteps[i] = i < sectorCount - 1 ? (MatrixXs)(M[i].transpose()) * steps[i + 1] : steps[i + 1];
        if (i + 1 < sectorCount)
            U[i + 1] = conditionY(endSteps[i]);
        steps[i] = sectors[i]->propagate(E, endSteps[i], sectors[i]->max, sectors[i]->min, true);
    }
    const Y<Scalar, Dynamic> matchRight = steps[matchIndex + 1];

    U[0] = MatrixXs::Identity(basisSize, basisSize);
    for (int i = 0; i <= matchIndex; ++i) {
        endSteps[i] = sectors[i]->propagate(E, steps[i], sectors[i]->min, sectors[i]->max, true);
        steps[i + 1] = M[i] * endSteps[i];
        U[i + 1] = conditionY(steps[i + 1]);
    }
    const Y<Scalar, Dynamic> matchLeft = steps[matchIndex + 1];
    MatrixXs kernel = getKernel<Scalar>(matchLeft, matchRight, 1e-5);

    vector<Y<Scalar, Dynamic>> elements;
    if (kernel.cols() > 0) {
        elements.resize(sectorCount + 1);
        MatrixXs left = matchLeft.getY(0).colPivHouseholderQr().solve(kernel);
        MatrixXs right = matchRight.getY(0).colPivHouseholderQr().solve(kernel);

        Y<Scalar, Dynamic> elementMatchRight = matchRight * right;

        MatrixXs normalizer = MatrixXs::Zero(left.cols(), left.cols());
        for (int i = matchIndex + 1; i >= 0; --i) {
            elements[i] = steps[i] * left;
            if (i <= matchIndex) {
                cec_cce<Scalar, true>(normalizer, endSteps[i] * left);
                cec_cce<Scalar, false>(normalizer, elements[i]);
            }
            if (i > 0) {
                U[i].template triangularView<Upper>().
                        template solveInPlace<OnTheLeft>(left);
            }
        }

        for (int i = matchIndex + 1; i < sectorCount; ++i) {
            elements[static_cast<size_t>(i + 1)] = endSteps[i] * right;
            cec_cce<Scalar, true>(normalizer, elements[i + 1]);
            cec_cce<Scalar, false>(normalizer, i > matchIndex + 1 ? steps[i] * right : elementMatchRight);
            if (i + 1 < sectorCount) {
                U[i + 1].template triangularView<Upper>().
                        template solveInPlace<OnTheLeft>(right);
            }
        }

        LLT<Ref<MatrixXs>> llt(normalizer); // In place

        for (Index i = 0; i <= sectorCount; ++i) {
            Y<Scalar, Dynamic> &element = elements[static_cast<size_t>(i)];
            llt.matrixL().transpose().template solveInPlace<OnTheRight>(element.data);
        }
    }
    delete[] steps;
    delete[] endSteps;
    delete[] U;
    return elements;
}

#include "../util/instantiate.h"
