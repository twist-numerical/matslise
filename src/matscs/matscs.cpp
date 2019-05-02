#include <functional>
#include "../matscs.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;

int find_sector(const Matscs *ms, double point) {
    int a = 0, b = ms->sectorCount, c;
    while (!ms->sectors[c = a + (b - a) / 2]->contains(point)) {
        if (c == a)
            return -1;
        if (point < ms->sectors[c]->min)
            b = c;
        else
            a = c;
    }
    return c;
}

template<int r>
Y<Dynamic, r>
Matscs::propagate(double E, const Y<Dynamic, r> &_y, double a, double b, bool use_h) const {
    if (!contains(a) || !contains(b))
        throw runtime_error("Matslise::propagate(): a and b should be in the interval");
    Y<Dynamic, r> y = _y;
    int sectorIndex = find_sector(this, a);
    int direction = a < b ? 1 : -1;
    Sector *sector;
    do {
        sector = sectors[sectorIndex];
        y = sector->propagate(E, y, a, b, use_h);
        sectorIndex += direction;
    } while (!sector->contains(b));
    return y;
}

MatrixXd Matscs::propagatePsi(double E, const MatrixXd &_psi, double a, double b) const {
    MatrixXd psi = _psi;
    if (a < b) {
        for (int i = 0; i < sectorCount; ++i) {
            Sector *sector = sectors[i];
            if (sector->max > a) {
                if (sector->min < a) // first
                    psi = sector->propagatePsi(E, psi, sector->min - a);

                if (sector->max >= b) { // last
                    psi = sector->propagatePsi(E, psi, b - sector->min);
                    break;
                }

                psi = sector->propagatePsi(E, psi, sector->max - sector->min);
            }
        }
    } else {
        for (int i = sectorCount - 1; i >= 0; --i) {
            Sector *sector = sectors[i];
            if (sector->min < a) {
                if (sector->max > a) // first
                    psi = sector->propagatePsi(E, psi, sector->min - a);
                else
                    psi = sector->propagatePsi(E, psi, sector->min - sector->max);

                if (sector->min <= b) { // last
                    psi = sector->propagatePsi(E, psi, b - sector->min);
                    break;
                }

            }
        }
    }
    return psi;
}

template Y<-1, -1>
Matscs::propagate<-1>(double E, const Y<-1, -1> &y, double a, double b, bool use_h) const;

template Y<-1, 1>
Matscs::propagate<1>(double E, const Y<-1, 1> &y, double a, double b, bool use_h) const;

Matscs::~Matscs() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
}

vector<Y<Dynamic>> *Matscs::computeEigenfunction(double E, vector<double> &x) const {
    sort(x.begin(), x.end());
    vector<Y<Dynamic>> *ys = new vector<Y<Dynamic>>();

    auto iterator = x.begin();

    while (iterator != x.end() && *iterator < xmin - EPS)
        iterator = x.erase(iterator);

    Sector *sector;
    Y<Dynamic> y(n);
    for (int i = 0; iterator != x.end(); ++iterator) {
        while ((sector = sectors[i])->max < *iterator) {
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