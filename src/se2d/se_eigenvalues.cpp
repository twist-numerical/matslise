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
vector<Scalar> SE2D<Scalar>::findEigenvalues(const Scalar &Emin, const Scalar &Emax) const {
    int startDepth = 3;
    int depth = 8;
    Scalar linear = 0.001;

    Scalar length = Emax - Emin;
    for (int i = 0; i < startDepth; ++i)
        length *= .5;

    set<Scalar> todo;
    for (Scalar a = Emin + length / 2; a < Emax - length / 4; a += length) {
        todo.insert(a);
    }

    for (int i = 0; i < depth; ++i) {
        set<Scalar> next;
        for (const Scalar &E : todo) {
            for (const auto &error : calculateErrors(E)) {
                Scalar d = error.first / error.second;
                if (abs(d) < length) {
                    next.insert(E - d);
                }
            }
        }

        todo.clear();
        auto from = next.begin();
        for (Scalar v = Emin; from != next.end(); v += length) {
            auto to = from;
            while (to != next.end() && *to < v) ++to;
            if (from != to) {
                todo.insert(mean<Scalar>(from, to));
                from = to;
            }
        }
        length *= .5;
    }

    vector<Scalar> values;
    for (const Scalar &guess : todo) {
        const Scalar E = findEigenvalue(guess);
        if (!isnan(E)) {
            for (const auto &error : calculateErrors(E)) {
                const Scalar d = error.first / error.second;
                if (error.first < 1 && abs(d) < linear) {
                    values.push_back(E - d);
                }
            }
        }
    }

    return values;
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