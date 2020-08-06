#include <iostream>
#include <map>
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
                               const matslise::Rectangle<2, Scalar> &domain,
                               const Options2<Scalar> &_options):
        AbstractMatslise2D<Scalar>(potential, domain), N(_options._N), options(_options) {
    dirichletBoundary = Y<Scalar, Eigen::Dynamic>::Dirichlet(N);
    auto sectorsBuild = options._builder(this, domain.min, domain.max);
    sectors = std::move(sectorsBuild.sectors);
    matchIndex = sectorsBuild.matchIndex;
    sectorCount = sectors.size();
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
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
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
pair<Scalar, Scalar>
Matslise2D<Scalar>::matchingError(const Y<Scalar, Eigen::Dynamic> &yLeft, const Scalar &E, bool use_h) const {
    vector<pair<Scalar, Scalar>> errors = matchingErrors(yLeft, E, use_h);
    return *min_element(errors.begin(), errors.end(), [](
            const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) -> bool {
        if (abs(a.first) > 100 || abs(b.first) > 100)
            return abs(a.first) < abs(b.first);
        return abs(a.first / a.second) < abs(b.first / b.second);
    });
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
Scalar Matslise2D<Scalar>::eigenvalue(const Y<Scalar, Dynamic> &left, const Scalar &_E, bool use_h) const {
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
        return NAN;
    return E;
}

template<typename Scalar>
Scalar Matslise2D<Scalar>::eigenvalueError(const Y<Scalar, Dynamic> &left, const Scalar &E) const {
    return abs(E - eigenvalue(left, E, false));
}

template<typename Scalar, typename I>
Scalar mean(I start, I end) {
    Scalar sum = 0;
    int count = 0;
    for (I i = start; i != end; ++i) {
        sum += *i;
        ++count;
    }
    return sum / count;
}

template<typename Scalar>
inline bool set_contains(const set<Scalar> &values, const Scalar &guess, const Scalar &length) {
    auto it = values.lower_bound(guess - length);
    return it != values.end() && *it <= guess + length;
}

template<typename Scalar>
vector<Scalar> Matslise2D<Scalar>::eigenvalues(
        const Y<Scalar, Dynamic> &left, const Scalar &Emin, const Scalar &Emax) const {
    typedef pair<Scalar, Scalar> Interval;
    const int initialSteps = 16;

    const Scalar linear = 0.001;
    const Scalar minLength = 0.001;

    function<bool(Interval, Interval)> comp = [](Interval a, Interval b) -> bool {
        return a.first > b.first;
    };

    vector<Interval> intervals;
    Scalar length = (Emax - Emin) / initialSteps;
    {
        Scalar a = Emin + length / 2;
        while (a < Emax - length / 4) {
            intervals.push_back({length / 2, a});
            push_heap(intervals.begin(), intervals.end(), comp);
            a += length;
        }
    }

    set<Scalar> eigenvalues;
    Scalar guess;
    while (!intervals.empty()) {
        tie(length, guess) = intervals.front();
        pop_heap(intervals.begin(), intervals.end(), comp);
        intervals.pop_back();

        if (!set_contains(eigenvalues, guess, length)) {
            for (const pair<Scalar, Scalar> &error : matchingErrors(left, guess)) {
                Scalar d = error.first / error.second;
                if (abs(d) < minLength) {
                    Scalar E = eigenvalue(left, guess - d);
                    if (!isnan(E) && !set_contains(eigenvalues, E, minLength))
                        eigenvalues.insert(E);
                } else if (abs(d) < 2 * length) {
                    intervals.push_back({min(abs(d), .5 * length), guess - d});
                    push_heap(intervals.begin(), intervals.end(), comp);
                }
            }
        }
    }

    vector<Scalar> result;
    for (const Scalar &E : eigenvalues) {
        long valueCount = static_cast<long>(result.size());
        for (const auto &error : matchingErrors(left, E)) {
            const Scalar d = error.first / error.second;
            if (error.first < 1 && abs(d) < linear) {
                result.push_back(E - d);
            }
        }
        if (valueCount != static_cast<long>(result.size()))
            sort(result.begin() + valueCount, result.end());
    }
    return result;
}

template<typename Scalar>
inline bool
is_first_eigenvalue(const Matslise2D<Scalar> &se2d, const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) {
    vector<Eigenfunction2D<Scalar>> eigenfunctions = se2d.eigenfunction(
            left, E);
    if (eigenfunctions.size() != 1)
        return false;
    Array<Scalar, Dynamic, Dynamic> values = eigenfunctions[0](
            Matslise2D<Scalar>::ArrayXs::LinSpaced(101, se2d.domain.sub.min, se2d.domain.sub.max),
            Matslise2D<Scalar>::ArrayXs::LinSpaced(101, se2d.domain.min, se2d.domain.max));
    return values.minCoeff() * values.maxCoeff() > -1e-2;
}

template<typename Scalar>
constexpr Scalar fundamentalGap(const Rectangle<2, Scalar> &domain) {
    // https://arxiv.org/abs/1006.1686
    Scalar step = constants<Scalar>::PI / domain.diameter();
    step *= step;
    step *= 3;
    return step;
}

template<typename Scalar>
vector<Scalar> Matslise2D<Scalar>::firstEigenvalues(const Y<Scalar, Dynamic> &left, int n) const {
    Scalar E0 = firstEigenvalue(left);
    const Scalar step = max(Scalar(1), fundamentalGap(domain));
    vector<Scalar> values{E0};
    Scalar start = E0;
    while (values.size() < (unsigned long) n) {
        for (auto E : eigenvalues(left, start, start + step)) {
            auto lowerBound = lower_bound(values.begin(), values.end(), E);
            if (lowerBound == values.end()) {
                if (abs(lowerBound[-1] - E) > 1e-5)
                    values.emplace_back(E);
            } else if (abs(lowerBound[0] - E) > 1e-5) {
                if (lowerBound - 1 == values.begin() || abs(lowerBound[-1] - E) > 1e-5)
                    values.insert(lowerBound, E);
            }
        }
        start += step;
    }
    while (values.size() > (unsigned long) n)
        values.pop_back();
    return values;
}

template<typename Scalar>
int numberPositive(const vector<pair<Scalar, Scalar>> &errors) {
    int n = 0;
    for (auto &error : errors)
        if (error.first > 0)
            ++n;
    return n;
}

template<typename Scalar>
Scalar Matslise2D<Scalar>::firstEigenvalue(const Y<Scalar, Eigen::Dynamic> &left) const {
    Scalar minimal = estimatePotentialMinimum();
    Scalar lower = minimal;
    cout << "Potential minimum: " << lower << endl;
    Scalar upper = eigenvalue(left, lower);
    if (is_first_eigenvalue(*this, left, upper))
        return upper;

    {
        int steps = 0;
        Scalar stepSize = 1;
        while (numberPositive(matchingErrors(left, lower)) != N) {
            lower -= (stepSize *= 2);
            if (++steps > 10) {
                throw runtime_error("Could not find the first eigenvalue. A useful lower bound is not found.");
            }
        }
    }
    {
        int steps = 0;
        Scalar stepSize = 2;
        while (numberPositive(matchingErrors(left, upper)) >= N) {
            upper += (stepSize *= 2);
            if (++steps > 10) {
                throw runtime_error("Could not find the first eigenvalue. A useful upper bound is not found.");
            }
        }
    }

    while (upper - lower > 1e-3) {
        Scalar mid = (lower + upper) / 2;
        if (numberPositive(matchingErrors(left, mid)) == N)
            lower = mid;
        else
            upper = mid;
    }
    Scalar first = eigenvalue(left, lower);
    cout << "Guess for first: " << first << endl;
    if (is_first_eigenvalue(*this, left, first))
        return first;

    first -= fundamentalGap(domain);
    cout << "Wrong! Going for bruteforce between " << minimal << " and " << first << endl;
    Index n = 60;
    Scalar h = (first - minimal) / n;
    Scalar guess = minimal;
    for (Index i = 1; i < n; ++i) {
        guess += h;
        if (numberPositive(matchingErrors(left, guess)) != N) {
            guess = eigenvalue(left, guess);
            if (is_first_eigenvalue(*this, left, guess))
                return guess;
            break;
        }
    }

    throw runtime_error("Could not find the first eigenvalue.");
}

template<typename Scalar>
vector<Scalar> Matslise2D<Scalar>::eigenvaluesByIndex(const Y<Scalar, Dynamic> &left, int Imin, int Imax) const {
    vector<Scalar> eigenvalues = firstEigenvalues(left, Imax);
    eigenvalues.erase(eigenvalues.begin(), eigenvalues.begin() + Imin);
    return eigenvalues;
}

#include "../util/instantiate.h"
