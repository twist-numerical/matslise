#ifndef MATSLISE_SECTORBUILDER_H
#define MATSLISE_SECTORBUILDER_H

#include "constants.h"
#include <stdexcept>

namespace matslise {
    template<typename Problem>
    struct SectorBuilderReturn {
        std::vector<typename Problem::Sector *> sectors;
        // The index of the sector left of the matching point
        int matchIndex;
    };

    template<typename Problem, typename Scalar=typename Problem::Scalar>
    using SectorBuilder=std::function<SectorBuilderReturn<Problem>(
            const Problem *, const Scalar &min, const Scalar &max)>;

    namespace sector_builder {
        template<typename Problem>
        SectorBuilder<Problem> uniform(int sectorCount) {
            typedef typename Problem::Scalar Scalar;
            typedef typename Problem::Sector Sector;
            return [sectorCount](const Problem *problem, const Scalar &min, const Scalar &max)
                    -> SectorBuilderReturn<Problem> {
                Scalar h = (max - min) / sectorCount;

                if (sectorCount == 1) {
                    return {
                            {new Sector(problem, min, max, false)},
                            0
                    };
                }

                std::vector<Sector *> sectors;
                sectors.resize(sectorCount);

                sectors[0] = new Sector(problem, min, min + h, false);
                sectors[sectorCount - 1] = new Sector(problem, min + (sectorCount - 1) * h, max, true);
                int i = 0, j = sectorCount - 1;
                while (i + 1 != j) {
                    if (Sector::compare(*sectors[j], *sectors[i])) {
                        ++i;
                        sectors[i] = new Sector(problem, min + i * h, min + (i + 1) * h, false);
                    } else {
                        --j;
                        sectors[j] = new Sector(problem, min + j * h, min + (j + 1) * h, true);
                    }
                }
                int matchIndex = i;
                return {std::move(sectors), matchIndex};
            };
        }

        template<typename Problem>
        typename Problem::Sector *automaticNextSector(
                const Problem *ms, const typename Problem::Scalar &_h,
                const typename Problem::Scalar &left, const typename Problem::Scalar &right,
                const typename Problem::Scalar &tolerance, bool forward) {
            typedef typename Problem::Scalar Scalar;
            typedef typename Problem::Sector Sector;
            Scalar h = _h;
            Scalar xmin = forward ? left : right - h;
            Scalar xmax = forward ? left + h : right;
            if (forward && xmax > right) {
                xmax = right;
            } else if (!forward && xmin < left) {
                xmin = left;
            }
            h = xmax - xmin;
            Sector *s = nullptr;
            Scalar error = 1;
            int steps = 0;
            do {
                if (s != nullptr) {
                    ++steps;
                    h *= std::max(Scalar(.1), pow(tolerance / error, 1. / (Problem::order - 1)));
                    if (forward) {
                        xmax = xmin + h;
                    } else {
                        xmin = xmax - h;
                    }
                    delete s;
                }
                s = new Sector(ms, xmin, xmax, !forward);
                error = s->error();
                // cout << "(h: " << h << ", error: " << error << ") ";
            } while (error > tolerance && steps < 10 && h > 1e-5);
            if (steps == 0) {
                while (error < tolerance / 2 && steps < 10 && h != right - left) {
                    ++steps;
                    h *= pow(tolerance / error, 1. / (Problem::order - 1));
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
                    auto *newSector = new Sector(ms, xmin, xmax, !forward);
                    error = newSector->error();
                    // cout << "(h: " << h << ", error: " << error << ") ";
                    if (error > tolerance) {
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

        template<typename Problem>
        SectorBuilder<Problem> automatic(const typename Problem::Scalar &tolerance) {
            typedef typename Problem::Scalar Scalar;
            typedef typename Problem::Sector Sector;
            return [tolerance](const Problem *problem, const Scalar &min, const Scalar &max)
                    -> SectorBuilderReturn<Problem> {
                std::vector<typename Problem::Sector *> forward;
                std::vector<typename Problem::Sector *> backward;
                // It shouldn't be the true middle
                //  so just dividing by something ~= 2
                typename Problem::Scalar mid = (max + min) / 2.11803398875;
                typename Problem::Scalar h = mid - min;
                forward.push_back(automaticNextSector(problem, h, min, mid, tolerance, true));
                backward.push_back(automaticNextSector(problem, h, mid, max, tolerance, false));


                while (forward.back()->max != backward.back()->min) {
                    if (Sector::compare(*backward.back(), *forward.back()))
                        forward.push_back(automaticNextSector(
                                problem, forward.back()->max - forward.back()->min, forward.back()->max,
                                backward.back()->min, tolerance, true));
                    else
                        backward.push_back(automaticNextSector(
                                problem, backward.back()->max - backward.back()->min, forward.back()->max,
                                backward.back()->min, tolerance, false));
                }

                int matchIndex = forward.size() - 1;
                forward.insert(forward.end(), backward.rbegin(), backward.rend());
                return {std::move(forward), matchIndex};
            };
        }


    }
}

#endif //MATSLISE_SECTORBUILDER_H
