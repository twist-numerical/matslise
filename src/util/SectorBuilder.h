#ifndef MATSLISE_SECTORBUILDER_H
#define MATSLISE_SECTORBUILDER_H

#include "SectorBuilder_header.h"
#include "../matslise.h"
#include "../matscs.h"
#include "../se2d.h"
#include "constants.h"
#include <stdexcept>

template<typename Problem>
inline bool matslise::sectorbuilder::compareSectors(
        const typename Problem::Sector &a, const typename Problem::Sector &b) {
    if (std::is_same<Problem, matslise::Matslise<typename Problem::Scalar>>::value) {
        return ((const typename Matslise<typename Problem::Scalar>::Sector &) a).vs[0] <
               ((const typename Matslise<typename Problem::Scalar>::Sector &) b).vs[0];
    } else if (std::is_same<Problem, matslise::Matscs<typename Problem::Scalar>>::value) {
        return (((const typename Matscs<typename Problem::Scalar>::Sector &) a).vs[0].diagonal() -
                ((const typename Matscs<typename Problem::Scalar>::Sector &) b).vs[0].diagonal()).sum() < 0;
    } else if (std::is_same<Problem, matslise::Matslise2D<typename Problem::Scalar>>::value) {
        return ((const typename Matslise2D<typename Problem::Scalar>::Sector &) a).vbar.minCoeff() <
               ((const typename Matslise2D<typename Problem::Scalar>::Sector &) b).vbar.minCoeff();
    }
    throw std::invalid_argument("Not supported");
}

template<typename Problem>
void matslise::sectorbuilder::Uniform<Problem>::build(
        Problem *ms, typename Problem::Scalar min, typename Problem::Scalar max) const {
    typename Problem::Scalar h = (max - min) / sectorCount;

    ms->sectors.resize(sectorCount);
    if (sectorCount == 1) {
        ms->sectors[0] = new typename Problem::Sector(ms, min, max, false);
        return;
    }

    ms->sectors[0] = new typename Problem::Sector(ms, min, min + h, false);
    ms->sectors[sectorCount - 1] = new typename Problem::Sector(ms, min + (sectorCount - 1) * h, max, true);
    int i = 0, j = sectorCount - 1;
    while (i + 1 != j) {
        if (compareSectors<Problem>(*ms->sectors[j], *ms->sectors[i])) {
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
void matslise::sectorbuilder::Auto<Problem>::build(
        Problem *ms, typename Problem::Scalar min, typename Problem::Scalar max) const {
    std::vector<typename Problem::Sector *> forward;
    std::vector<typename Problem::Sector *> backward;
    // It shouldn't be the true middle
    //  so just dividing by something ~= 2
    typename Problem::Scalar mid = (max + min) / 2.11803398875;
    typename Problem::Scalar h = mid - min;
    forward.push_back(nextSector<true>(ms, h, min, mid));
    backward.push_back(nextSector<false>(ms, h, mid, max));


    while (forward.back()->max != backward.back()->min) {
        if (compareSectors<Problem>(*backward.back(), *forward.back()))
            forward.push_back(nextSector<true>(ms, forward.back()->max - forward.back()->min,
                                               forward.back()->max, backward.back()->min));
        else
            backward.push_back(nextSector<false>(ms, backward.back()->max - backward.back()->min,
                                                 forward.back()->max, backward.back()->min));
    }

    ms->match = forward.back()->max;
    ms->sectors.reserve(forward.size() + backward.size());
    for (typename Problem::Sector *const &s : forward)
        ms->sectors.push_back(s);
    for (auto j = backward.rbegin(); j != backward.rend(); ++j)
        ms->sectors.push_back(*j);


    /*  for (int i = 0; i < ms->sectorCount; ++i)
          cout << "h: " << (ms->sectors[i]->xmax - ms->sectors[i]->xmin) << " (" << ms->sectors[i]->xmin << ", "
               << ms->sectors[i]->xmax << ") error: " << ms->sectors[i]->calculateError() << endl;
      cout << "match: " << ms->match << "\n" << endl;*/
}


template<typename Problem>
template<bool forward>
typename Problem::Sector *matslise::sectorbuilder::Auto<Problem>::nextSector(
        const Problem *ms, const typename Problem::Scalar &_h, const typename Problem::Scalar &left,
        const typename Problem::Scalar &right) const {
    typename Problem::Scalar h = _h;
    typename Problem::Scalar xmin = forward ? left : right - h;
    typename Problem::Scalar xmax = forward ? left + h : right;
    if (forward && xmax > right) {
        xmax = right;
    } else if (!forward && xmin < left) {
        xmin = left;
    }
    h = xmax - xmin;
    typename Problem::Sector *s = nullptr;
    typename Problem::Scalar error = 1;
    int steps = 0;
    do {
        if (s != nullptr) {
            ++steps;
            h *= std::max((typename Problem::Scalar) .1, pow(tol / error, 1. / (Problem::order - 1)));
            if (forward) {
                xmax = xmin + h;
            } else {
                xmin = xmax - h;
            }
            delete s;
        }
        s = new typename Problem::Sector(ms, xmin, xmax, !forward);
        error = s->error();
        // cout << "(h: " << h << ", error: " << error << ") ";
    } while (error > tol && steps < 10 && h > 1e-5);
    if (steps == 0) {
        while (error < tol / 2 && steps < 10 && h != right - left) {
            ++steps;
            h *= pow(tol / error, 1. / (Problem::order - 1));
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
            auto *newSector = new typename Problem::Sector(ms, xmin, xmax, !forward);
            error = newSector->error();
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

#endif //MATSLISE_SECTORBUILDER_H
