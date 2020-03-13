#include "../matslise.h"
#include <set>

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;

template<typename Scalar>
Scalar SE2D<Scalar>::findEigenvalue(
        const Scalar &_E, const Scalar &tolerance, int maxIterations, const Scalar &minTolerance) const {
    Scalar E = _E;
    Scalar error, derror;
    int i = 0;
    do {
        tie(error, derror) = calculateError(E, &NEWTON_RAPHSON_SORTER<Scalar>);
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
vector<Scalar> SE2D<Scalar>::findEigenvalues(const Scalar &Emin, const Scalar &Emax, const int &initialSteps) const {
    typedef pair<Scalar, Scalar> Interval;

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
            for (const pair<Scalar, Scalar> &error : calculateErrors(guess)) {
                Scalar d = error.first / error.second;
                if (abs(d) < minLength) {
                    Scalar E = findEigenvalue(guess - d);
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
        for (const auto &error : calculateErrors(E)) {
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
Scalar SE2D<Scalar>::findFirstEigenvalue() const {
    // TODO: still a WIP
    Scalar E = 0;
    bool changed;
    Scalar err, derr;
    do {
        E = findEigenvalue(0);
        cout << E << endl;
        changed = false;
        for (const pair<Scalar, Scalar> &errorPair: calculateErrors(E)) {
            tie(err, derr) = errorPair;
            Scalar newE = E - err / derr;
            cout << "       " << newE << endl;
            if (newE < E) {
                cout << E << " >> " << newE << endl;
                E = newE;
                changed = abs(E - newE) > 1e-9;
            }

        }
    } while (changed);
    return E;
}

template<typename Scalar>
vector<Scalar> SE2D<Scalar>::computeEigenvaluesByIndex(int Imin, int Imax) const {
    // TODO
    (void) Imin;
    (void) Imax;
    vector<Scalar> values;
    return values;
}

#include "../util/instantiate.h"