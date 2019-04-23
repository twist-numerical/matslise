#include <functional>
#include "../matscs.h"

#define EPS (1.e-12)

using namespace std;
using namespace matslise;

template<int r>
Y<Dynamic, r>
Matscs::propagate(double E, const Y<Dynamic, r> &_y, double a, double b, bool use_h) const {
    Y<Dynamic, r> y = _y;
    if (a < b) {
        for (int i = 0; i < sectorCount; ++i) {
            Sector *sector = sectors[i];
            if (sector->xmax > a) {
                if (sector->xmin < a) // first
                    y = sector->propagate(E, y, sector->xmin - a, use_h);

                if (sector->xmax >= b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, use_h);
                    break;
                }

                y = sector->propagate(E, y, sector->xmax - sector->xmin, use_h);
            }
        }
    } else {
        for (int i = sectorCount - 1; i >= 0; --i) {
            Sector *sector = sectors[i];
            if (sector->xmin < a) {
                if (sector->xmax > a) // first
                    y = sector->propagate(E, y, sector->xmin - a, use_h);
                else
                    y = sector->propagate(E, y, sector->xmin - sector->xmax, use_h);

                if (sector->xmin <= b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, use_h);
                    break;
                }

            }
        }
    }
    return y;
}

MatrixXd Matscs::propagatePsi(double E, const MatrixXd &_psi, double a, double b) const {
    MatrixXd psi = _psi;
    if (a < b) {
        for (int i = 0; i < sectorCount; ++i) {
            Sector *sector = sectors[i];
            if (sector->xmax > a) {
                if (sector->xmin < a) // first
                    psi = sector->propagatePsi(E, psi, sector->xmin - a);

                if (sector->xmax >= b) { // last
                    psi = sector->propagatePsi(E, psi, b - sector->xmin);
                    break;
                }

                psi = sector->propagatePsi(E, psi, sector->xmax - sector->xmin);
            }
        }
    } else {
        for (int i = sectorCount - 1; i >= 0; --i) {
            Sector *sector = sectors[i];
            if (sector->xmin < a) {
                if (sector->xmax > a) // first
                    psi = sector->propagatePsi(E, psi, sector->xmin - a);
                else
                    psi = sector->propagatePsi(E, psi, sector->xmin - sector->xmax);

                if (sector->xmin <= b) { // last
                    psi = sector->propagatePsi(E, psi, b - sector->xmin);
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
        while ((sector = sectors[i])->xmax < *iterator) {
            y = sector->calculateT(E) * y;
            ++i;
            if (i >= sectorCount)
                goto allSectorsDone;
        }

        ys->push_back(sector->calculateT(E, *iterator - sector->xmin) * y);
    }

    allSectorsDone:
    while (iterator != x.end() && *iterator > xmax + EPS)
        iterator = x.erase(iterator);

    return ys;
}