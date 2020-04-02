#include <functional>
#include "../matslise.h"
#include "../util/find_sector.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;
using namespace Eigen;

template<typename Scalar>
template<int r>
Y<Scalar, Dynamic, r>
Matscs<Scalar>::propagateColumn(const Scalar &E, const Y<Scalar, Dynamic, r> &_y,
                                const Scalar &a, const Scalar &b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Scalar, Dynamic, r> y = _y;
    int sectorIndex = find_sector<Matscs<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        y = sector->propagateColumn(E, y, a, b, use_h);
        sectorIndex += direction;
    } while (!sector->contains(b));
    return y;
}

template<typename Scalar>
pair<Y<Scalar, Dynamic>, Scalar>
Matscs<Scalar>::propagate(const Scalar &E, const Y<Scalar, Dynamic> &_y,
                          const Scalar &a, const Scalar &b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matscs::propagate(): a and b should be in the interval");
    Y<Scalar, Dynamic> y = _y;
    int sectorIndex = find_sector<Matscs<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    Scalar dTheta;
    Scalar argdet = 0;
    do {
        sector = sectors[sectorIndex];
        tie(y, dTheta) = sector->propagate(E, y, a, b, use_h);
        argdet += dTheta;
        sectorIndex += direction;
    } while (!sector->contains(b));
    return {y, argdet};
}

template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic> Matscs<Scalar>::propagatePsi(
        const Scalar &E, const Matrix<Scalar, Dynamic, Dynamic> &_psi, const Scalar &a, const Scalar &b) const {
    Matrix<Scalar, Dynamic, Dynamic> psi = _psi;

    int sectorIndex = find_sector<Matscs<Scalar>>(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        psi = sector->propagatePsi(E, psi, a, b);
        sectorIndex += direction;
    } while (!sector->contains(b));
    return psi;
}

template<typename Scalar>
Matscs<Scalar>::~Matscs() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
}

template<typename Scalar>
vector<Y<Scalar, Dynamic>> *Matscs<Scalar>::eigenfunction(const Scalar &E, vector<Scalar> &x) const {
    sort(x.begin(), x.end());
    auto *ys = new vector<Y<Scalar, Dynamic>>();

    auto iterator = x.begin();

    while (iterator != x.end() && *iterator < xmin - EPS)
        iterator = x.erase(iterator);

    Y<Scalar, Dynamic> y(n);
    for (int i = 0; iterator != x.end(); ++iterator) {
        Sector *const &sector = sectors[i];
        while (*iterator > sector->max) {
            y = sector->calculateT(E) * y;
            ++i;
            if (i >= sectorCount)
                goto allSectorsDone;
        }

        ys->push_back(sector->calculateT(E, *iterator - sector->min) * y);
    }

    allSectorsDone:
    while (iterator != x.end() && *iterator > xmax + EPS)
        iterator = x.erase(iterator);

    return ys;
}

#define INSTANTIATE_PROPAGATE(Scalar, r)\
template Y<Scalar, Dynamic, r>\
Matscs<Scalar>::propagateColumn<r>(\
    const Scalar &E, const Y<Scalar, Dynamic, r> &y, const Scalar &a, const Scalar &b, bool use_h) const;

#define INSTANTIATE_MORE(Scalar)\
INSTANTIATE_PROPAGATE(Scalar, 1)\
INSTANTIATE_PROPAGATE(Scalar, -1)

#include "../util/instantiate.h"