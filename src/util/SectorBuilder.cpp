//
// Created by toon on 4/23/19.
//

#include "SectorBuilder.h"
#include "../matslise.h"
#include "../matscs.h"
#include "../se2d.h"

using namespace matslise;
using namespace matslise::sectorbuilder;
using namespace std;

template<typename Problem>
// a < b ?
bool compareSectors(typename Problem::Sector *a, typename Problem::Sector *b);

template<>
bool compareSectors<Matslise>(Matslise::Sector *a, Matslise::Sector *b) {
    return a->vs[0] < b->vs[0];
}

template<>
bool compareSectors<Matscs>(Matscs::Sector *a, Matscs::Sector *b) {
    return (a->vs[0].diagonal() - b->vs[0].diagonal()).sum() < 0;
}

template<>
bool compareSectors<SEnD<2>>(SEnD<2>::Sector *a, SEnD<2>::Sector *b) {
    return a->vbar.minCoeff() < b->vbar.minCoeff();
}

template<typename Problem>
void Uniform<Problem>::build(Problem *ms, double min, double max) const {
    ms->sectorCount = sectorCount;
    ms->sectors = new typename Problem::Sector *[sectorCount];
    double h = (max - min) / sectorCount;

    if (sectorCount == 1) {
        ms->sectors[0] = new typename Problem::Sector(ms, min, max, false);
        return;
    }

    ms->sectors[0] = new typename Problem::Sector(ms, min, min + h, false);
    ms->sectors[sectorCount - 1] = new typename Problem::Sector(ms, min + (sectorCount - 1) * h, max, true);
    int i = 0, j = sectorCount - 1;
    while (i + 1 != j) {
        if (compareSectors<Problem>(ms->sectors[j], ms->sectors[i])) {
            ++i;
            ms->sectors[i] = new typename Problem::Sector(ms, min + i * h, min + (i + 1) * h, false);
        } else {
            --j;
            ms->sectors[j] = new typename Problem::Sector(ms, min + j * h, min + (j + 1) * h, true);
        }
    }
    int matchIndex = i;

    ms->match = ms->sectors[matchIndex]->max;
}

template<typename Problem>
void Auto<Problem>::build(Problem *ms, double min, double max) const {
    vector<typename Problem::Sector *> forward;
    vector<typename Problem::Sector *> backward;
    double mid = (max + min) / 2;
    double h = mid - min;
    forward.push_back(nextSector<true>(ms, h, min, mid));
    backward.push_back(nextSector<false>(ms, h, mid, max));


    while (forward.back()->max != backward.back()->min) {
        if (compareSectors<Problem>(backward.back(), forward.back()))
            forward.push_back(nextSector<true>(ms, forward.back()->max - forward.back()->min,
                                               forward.back()->max, backward.back()->min));
        else
            backward.push_back(nextSector<false>(ms, backward.back()->max - backward.back()->min,
                                                 forward.back()->max, backward.back()->min));
    }

    ms->match = forward.back()->max;
    ms->sectorCount = (int) (forward.size() + backward.size());
    ms->sectors = new typename Problem::Sector *[ms->sectorCount];
    int i = 0;
    for (typename Problem::Sector *s : forward)
        ms->sectors[i++] = s;
    for (auto j = backward.rbegin(); j != backward.rend(); ++j)
        ms->sectors[i++] = *j;


    /*  for (int i = 0; i < ms->sectorCount; ++i)
          cout << "h: " << (ms->sectors[i]->xmax - ms->sectors[i]->xmin) << " (" << ms->sectors[i]->xmin << ", "
               << ms->sectors[i]->xmax << ") error: " << ms->sectors[i]->calculateError() << endl;
      cout << "match: " << ms->match << "\n" << endl;*/
}

template<typename Problem>
template<bool forward>
typename Problem::Sector *
Auto<Problem>::nextSector(Problem *ms, double h, double left, double right) const {
    double xmin = forward ? left : right - h;
    double xmax = forward ? left + h : right;
    if (forward && xmax > right) {
        xmax = right;
    } else if (!forward && xmin < left) {
        xmin = left;
    }
    h = xmax - xmin;
    typename Problem::Sector *s = nullptr;
    double error = 1;
    int steps = 0;
    do {
        if (s != nullptr) {
            ++steps;
            h *= max(.01, pow(tol / error, 1. / 6));
            if (forward) {
                xmax = xmin + h;
            } else {
                xmin = xmax - h;
            }
            delete s;
        }
        s = new typename Problem::Sector(ms, xmin, xmax, !forward);
        error = s->calculateError();
        // cout << "(h: " << h << ", error: " << error << ") ";
    } while (error > tol && steps < 10 && h > 1e-5);
    if (steps == 0) {
        while (error < tol / 2 && steps < 10 && h != right - left) {
            ++steps;
            h *= pow(tol / error, 1. / 12);
            if (forward) {
                xmax = xmin + h;
                if (xmax > right)
                    xmax = right;
            } else {
                xmin = xmax - h;
                if (xmin < left)
                    xmin = left;
            }
            h = xmax - xmin;
            typename Problem::Sector *newSector = new typename Problem::Sector(ms, xmin, xmax, !forward);
            error = newSector->calculateError();
            // cout << "(h: " << h << ", error: " << error << ") ";
            if (error > tol) {
                delete newSector;
                break;
            } else {
                delete s;
                s = newSector;
            }
        }
    }
    // cout << "-> " << h << endl;
    return s;
}

template void Uniform<Matslise>::build(Matslise *, double, double) const;

template void Auto<Matslise>::build(Matslise *, double, double) const;

template void Uniform<Matscs>::build(Matscs *, double, double) const;

template void Auto<Matscs>::build(Matscs *, double, double) const;

template void Uniform<SEnD<2>>::build(SEnD<2> *, double, double) const;

template void Auto<SEnD<2>>::build(SEnD<2> *, double, double) const;

