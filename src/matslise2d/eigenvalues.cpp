#include "../matslise.h"
#include <set>

using namespace Eigen;
using namespace matslise;
using namespace std;

template<typename Scalar>
Scalar Matslise2D<Scalar>::eigenvalue(const Y<Scalar, Dynamic> &left, const Scalar &_E) const {
    const Scalar tolerance = 1e-9;
    const Scalar minTolerance = 1e-5;
    const int maxIterations = 30;

    Scalar E = _E;
    Scalar error, derror;
    int i = 0;
    do {
        tie(error, derror) = matchingError(left, E);
        E -= error / derror;
        ++i;
    } while (i < maxIterations && abs(error) > tolerance);

    if (abs(error) > minTolerance)
        return NAN;
    return E;
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
    vector<typename Matslise2D<Scalar>::ArrayXXs> eigenfunctions = se2d.eigenfunction(
            left, E,
            Matslise2D<Scalar>::ArrayXs::LinSpaced(50, se2d.domain.sub.min, se2d.domain.sub.max),
            Matslise2D<Scalar>::ArrayXs::LinSpaced(50, se2d.domain.min, se2d.domain.max));
    if (eigenfunctions.size() != 1)
        return false;

    return eigenfunctions[0].minCoeff() * eigenfunctions[0].abs().maxCoeff() > -1e-5;
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
Scalar Matslise2D<Scalar>::firstEigenvalue(const Y<Scalar, Eigen::Dynamic> &left) const {
    Scalar potentialMinimum = sectors[0]->matslise->estimatePotentialMinimum();
    for (int i = 0; i < sectorCount; ++i)
        potentialMinimum = min(potentialMinimum, sectors[i]->matslise->estimatePotentialMinimum());

    const int maxDepth = 4;
    vector<pair<Scalar, int>> guesses{{potentialMinimum, 0}};
    function<bool(const pair<Scalar, int> &, const pair<Scalar, int> &)> comp
            = [](const pair<Scalar, int> &a, const pair<Scalar, int> &b) -> bool {
                return b.first < a.first;
            };

    int depth;
    Scalar guess;
    set<Scalar> seen;

    while (!guesses.empty()) {
        tie(guess, depth) = guesses.front();
        pop_heap(guesses.begin(), guesses.end(), comp);
        guesses.pop_back();

        if (depth >= maxDepth) {
            Scalar E = eigenvalue(left, guess);

            if (!set_contains(seen, E, 1e-5) && is_first_eigenvalue(*this, left, E))
                return E;
            seen.insert(E);
        } else {
            auto errors = matchingErrors(left, guess);
            sort(errors.begin(), errors.end(),
                 [](const pair<Scalar, Scalar> &a, const pair<Scalar, Scalar> &b) -> bool {
                     return a.first/a.second > b.first/b.second;
                 });
            for (auto it = errors.begin(); it != errors.end() && it != errors.begin() + 3; ++it) {
                Scalar newGuess = guess - it->first / it->second;
                guesses.emplace_back(newGuess, depth + 1);
                push_heap(guesses.begin(), guesses.end(), comp);
            }
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