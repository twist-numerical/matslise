//
// Created by toon on 4/23/19.
//

#include "SectorBuilder.h"
#include "../matslise.h"

using namespace matslise;
using namespace std;

template<>
void SectorBuilder<Matslise>::Uniform::build(Matslise *ms) const {
    ms->sectorCount = sectorCount;
    ms->sectors = new Matslise::Sector *[sectorCount];
    double h = (ms->xmax - ms->xmin) / sectorCount;

    for (int i = 0; i < sectorCount; ++i) {
        double a = ms->xmin + i * h;
        double b = ms->xmax - (sectorCount - i - 1) * h;
        ms->sectors[i] = new Matslise::Sector(ms, a, b);
    }

    int matchIndex = 0;
    for (int i = 1; i < sectorCount - 1; ++i) {
        if (ms->sectors[i]->vs[0] < ms->sectors[matchIndex]->vs[0])
            matchIndex = i;
    }
    ms->match = ms->sectors[matchIndex]->xmax;
}

template<>
void SectorBuilder<Matslise>::Auto::build(Matslise *ms) const {
    vector < Matslise::Sector * > forward;
    vector < Matslise::Sector * > backward;
    double mid = (ms->xmax + ms->xmin) / 2;
    double h = mid - ms->xmin;
    forward.push_back(nextSector<true>(ms, h, ms->xmin, mid));
    backward.push_back(nextSector<false>(ms, h, mid, ms->xmax));


    while (forward.back()->xmax != backward.back()->xmin) {
        if (forward.back()->vs[0] >= backward.back()->vs[0])
            forward.push_back(nextSector<true>(ms, forward.back()->xmax - forward.back()->xmin,
                                               forward.back()->xmax, backward.back()->xmin));
        else
            backward.push_back(nextSector<false>(ms, backward.back()->xmax - backward.back()->xmin,
                                                 forward.back()->xmax, backward.back()->xmin));
    }

    ms->match = forward.back()->xmax;
    ms->sectorCount = (int) (forward.size() + backward.size());
    ms->sectors = new Matslise::Sector *[ms->sectorCount];
    int i = 0;
    for (Matslise::Sector *s : forward)
        ms->sectors[i++] = s;
    for (auto j = backward.rbegin(); j != backward.rend(); ++j)
        ms->sectors[i++] = *j;
/*

    for (int i = 0; i < ms->sectorCount; ++i)
        cout << "h: " << (ms->sectors[i]->xmax - ms->sectors[i]->xmin) << " (" << ms->sectors[i]->xmin << ", "
             << ms->sectors[i]->xmax << ") error: " << ms->sectors[i]->calculateError() << endl;
    cout << "match: " << ms->match << "\n" << endl;*/
}

template<>
template<bool forward>
Matslise::Sector *SectorBuilder<Matslise>::Auto::nextSector(Matslise *ms, double h, double left, double right) const {
    double xmin = forward ? left : right - h;
    double xmax = forward ? left + h : right;
    if (forward && xmax > right) {
        xmax = right;
    } else if (!forward && xmin < left) {
        xmin = left;
    }
    h = xmax - xmin;
    Matslise::Sector *s = nullptr;
    double error = 1;
    int steps = 0;
    do {
        if (s != nullptr) {
            ++steps;
            h *= pow(tol / error, 1. / 6);
            if (forward) {
                xmax = xmin + h;
            } else {
                xmin = xmax - h;
            }
            delete s;
        }
        s = new Matslise::Sector(ms, xmin, xmax);
        error = s->calculateError();
        // cout << "(h: " << h << ", error: " << error << ") ";
    } while (error > tol && steps < 5 && h > 1e-5);
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
            Matslise::Sector *newSector = new Matslise::Sector(ms, xmin, xmax);
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
