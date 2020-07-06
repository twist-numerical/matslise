#include "./sectorbuilder.h"
#include "../matslise.h"

using namespace matslise;

// Determine if Problem::Sector::refine is defined
template<typename Scalar, typename Problem>
typename std::integral_constant<bool, std::is_same<decltype(
std::declval<typename Problem::Sector>().refine(
        std::declval<Problem *>(), std::declval<Scalar>(), std::declval<Scalar>(), std::declval<bool>())
), typename Problem::Sector *>::value> __sectorHasRefine(int);

template<typename Scalar, typename Problem>
std::false_type __sectorHasRefine(...);

template<typename Scalar, typename Problem>
constexpr bool sectorHasRefine = decltype(__sectorHasRefine<Scalar, Problem>(0))::value;

template<typename Scalar, typename Problem>
typename std::enable_if<sectorHasRefine<Scalar, Problem>, typename Problem::Sector *>::type
refineSector(const typename Problem::Sector &sector, const Problem *problem, const Scalar &min,
             const Scalar &max, Direction direction) {
    return sector.refine(problem, min, max, direction);
}


template<typename Scalar, typename Problem>
static typename std::enable_if<!sectorHasRefine<Scalar, Problem>, typename Problem::Sector *>::type
refineSector(const typename Problem::Sector &, const Problem *problem, const Scalar &min, const Scalar &max,
             Direction direction) {
    return new typename Problem::Sector(problem, min, max, direction);
}

template<typename Problem>
SectorBuilder<Problem> sector_builder::uniform(int sectorCount) {
    typedef typename Problem::Scalar Scalar;
    typedef typename Problem::Sector Sector;
    return [sectorCount](const Problem *problem, const Scalar &min, const Scalar &max)
            -> SectorBuilderReturn<Problem> {
        Scalar h = (max - min) / sectorCount;

        if (sectorCount == 1) {
            return {{new Sector(problem, min, max, forward)}, 0};
        }

        std::vector<Sector *> sectors;
        sectors.resize(sectorCount);
#ifdef MATSLISE_parallel
#pragma omp parallel for shared(min, h)
        for (int i = 0; i < sectorCount; ++i) {
            sectors[i] = new Sector(problem, min + i * h, max - (sectorCount - i - 1) * h, none);
        }
        int matchIndex = 0;
        for (int i = 1; i < sectorCount; ++i) {
            if (Sector::compare(*sectors[i], *sectors[matchIndex]))
                matchIndex = i-1;
        }

#pragma omp parallel for
        for (int i = 0; i < sectorCount; ++i) {
            sectors[i]->setDirection(i <= matchIndex ? forward : backward);
        }

#else
        sectors[0] = new Sector(problem, min, min + h, forward);
        sectors[sectorCount - 1] = new Sector(problem, min + (sectorCount - 1) * h, max, backward);
        int i = 0, j = sectorCount - 1;
        while (i + 1 != j) {
            if (Sector::compare(*sectors[j], *sectors[i])) {
                ++i;
                sectors[i] = new Sector(problem, min + i * h, min + (i + 1) * h, forward);
            } else {
                --j;
                sectors[j] = new Sector(problem, min + j * h, min + (j + 1) * h, backward);
            }
        }
        int matchIndex = i;
#endif
        return {std::move(sectors), matchIndex};
    };
}

template<typename Problem>
typename Problem::Sector *automaticNextSector(
        const Problem *ms, typename Problem::Scalar &h,
        const typename Problem::Scalar &left, const typename Problem::Scalar &right,
        const typename Problem::Scalar &tolerance, Direction direction) {
    typedef typename Problem::Scalar Scalar;
    typedef typename Problem::Sector Sector;
    Scalar xmin = direction == forward ? left : right - h;
    Scalar xmax = direction == forward ? left + h : right;
    if (direction == forward && xmax > right) {
        xmax = right;
    } else if (direction == backward && xmin < left) {
        xmin = left;
    }
    h = xmax - xmin;
    int steps = 0;
    Sector *s = new Sector(ms, xmin, xmax, direction);
    Scalar error = s->error();
    while (error > tolerance && steps < 10 && h > 1e-3) {
        ++steps;
        h *= std::max(Scalar(.1), Scalar(pow(tolerance / error, Scalar(1. / (Problem::order - 1)))));
        if (direction == forward) {
            xmax = xmin + h;
        } else {
            xmin = xmax - h;
        }
        Sector *prevSector = s;
        s = refineSector<Scalar, Problem>(*prevSector, ms, xmin, xmax, direction);
        h = xmax - xmin;
        delete prevSector;
        error = s->error();
    }
    if constexpr(Sector::expensive) {
        if (error < tolerance / 2) {
            if (error <= 0) {
                h = right - left;
            } else {
                h *= pow(tolerance / error, 1. / (Problem::order - 1));
                if (h > right - left)
                    h = right - left;
            }
        }
    } else {
        if (steps == 0) {
            while (error < tolerance / 2 && steps < 10 && h != right - left) {
                ++steps;
                if (error <= 0) {
                    // rare edge cases
                    xmin = left;
                    xmax = right;
                } else {
                    h *= pow(tolerance / error, 1. / (Problem::order - 1));
                    if (direction == forward) {
                        xmax = xmin + h;
                        if (xmax > right)
                            xmax = right;
                    } else {
                        xmin = xmax - h;
                        if (xmin < left)
                            xmin = left;
                    }
                }
                h = xmax - xmin;
                auto *newSector = refineSector<Scalar, Problem>(
                        *s, ms, xmin, xmax, direction);
                error = newSector->error();

                if (error > tolerance) {
                    h = s->max - s->min;
                    delete newSector;
                    break;
                } else {
                    delete s;
                    s = newSector;
                }
            }
        }
    }

    return s;
}

template<typename Problem>
SectorBuilder<Problem> sector_builder::automatic(const typename Problem::Scalar &tolerance) {
    typedef typename Problem::Scalar Scalar;
    typedef typename Problem::Sector Sector;
    return [tolerance](const Problem *problem, const Scalar &min, const Scalar &max)
            -> SectorBuilderReturn<Problem> {
        std::vector<typename Problem::Sector *> forward;
        std::vector<typename Problem::Sector *> backward;
        // It shouldn't be the true middle
        Scalar mid = 0.4956864123 * max + 0.5043135877 * min;
        Scalar forwardH = .33 * (max - min);
        Scalar backwardH = forwardH;
        forward.push_back(automaticNextSector(problem, forwardH, min, mid, tolerance, Direction::forward));
        backward.push_back(automaticNextSector(problem, backwardH, mid, max, tolerance, Direction::backward));


        while (forward.back()->max != backward.back()->min) {
            if (Sector::compare(*backward.back(), *forward.back()))
                forward.push_back(automaticNextSector(
                        problem, forwardH, forward.back()->max,
                        backward.back()->min, tolerance, Direction::forward));
            else
                backward.push_back(automaticNextSector(
                        problem, backwardH, forward.back()->max,
                        backward.back()->min, tolerance, Direction::backward));
        }

        int matchIndex = forward.size() - 1;
        forward.insert(forward.end(), backward.rbegin(), backward.rend());
        return {std::move(forward), matchIndex};
    };
}

#define INSTANTIATE_SECTOR_BUILDER(Problem) \
template SectorBuilder<Problem> sector_builder::uniform<Problem>(int); \
template SectorBuilder<Problem> sector_builder::automatic<Problem>(const Problem::Scalar&);

#define INSTANTIATE_MORE(Scalar) \
INSTANTIATE_SECTOR_BUILDER(Matslise<Scalar>) \
INSTANTIATE_SECTOR_BUILDER(Matscs<Scalar>) \
INSTANTIATE_SECTOR_BUILDER(Matslise2D<Scalar>) \
INSTANTIATE_SECTOR_BUILDER(Matslise3D<Scalar>)

#include "./instantiate.h"
