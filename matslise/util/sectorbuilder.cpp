#include "../matslise.h"
#include "./sectorbuilder.h"
#include <algorithm>
#include <tuple>


using namespace matslise;
using namespace matslise::sector_builder;
using std::unique_ptr;
using std::vector;

template<typename Scalar, typename Problem>
static typename Problem::Sector *
refineSector(const typename Problem::Sector &sector, const Problem *problem, const Scalar &min, const Scalar &max,
             Direction direction) {
    if constexpr(Problem::refineSectors) {
        return sector.refine(problem, min, max, direction);
    } else {
        return new typename Problem::Sector(problem, min, max, direction);
    }
}

template<typename Problem, typename Scalar>
SectorBuilderReturn<Problem>
UniformSectorBuilder<Problem, Scalar>::operator()(const Problem *problem, const Scalar &min, const Scalar &max) const {
    typedef typename Problem::Sector Sector;
    Scalar h = (max - min) / sectorCount;

    SectorBuilderReturn<Problem> r;
    auto &sectors = r.sectors;
    sectors.resize(sectorCount);

    if (sectorCount == 1) {
        sectors[0].reset(new Sector(problem, min, max, forward));
        r.matchIndex = 0;
    } else {
        sectors[0].reset(new Sector(problem, min, min + h, forward));
        sectors[sectorCount - 1].reset(new Sector(problem, min + (sectorCount - 1) * h, max, backward));
        int i = 0, j = sectorCount - 1;
        while (i + 1 != j) {
            if (Sector::compare(*sectors[j], *sectors[i])) {
                ++i;
                sectors[i].reset(new Sector(problem, min + i * h, min + (i + 1) * h, forward));
            } else {
                --j;
                sectors[j].reset(new Sector(problem, min + j * h, min + (j + 1) * h, backward));
            }
        }
        r.matchIndex = i;
    }

    return r;
}

template<typename Problem>
typename Problem::Sector *automaticNextSector(
        const Problem *ms, typename Problem::Scalar &h,
        const typename Problem::Scalar &left, const typename Problem::Scalar &right,
        const typename Problem::Scalar &tolerance, Direction direction) {
    typedef typename Problem::Scalar Scalar;
    typedef typename Problem::Sector Sector;
    Scalar xmin = direction == Direction::forward ? left : right - h;
    Scalar xmax = direction == Direction::forward ? left + h : right;
    if (direction == Direction::forward && xmax > right) {
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
        h *= std::max(Scalar(.1), .99 * Scalar(pow(tolerance / error, Scalar(1. / (Problem::order - 1)))));
        if (direction == Direction::forward) {
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
                    if (direction == Direction::forward) {
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

template<typename Problem, typename Scalar>
SectorBuilderReturn<Problem>
AutomaticSectorBuilder<Problem, Scalar>::operator()(const Problem *problem, const Scalar &min,
                                                    const Scalar &max) const {
    typedef typename Problem::Sector Sector;

    std::vector<unique_ptr<Sector>> forward;
    std::vector<unique_ptr<Sector>> backward;
    // It shouldn't be the true middle
    Scalar mid = 0.4956864123 * max + 0.5043135877 * min;
    Scalar forwardH = .33 * (max - min);
    Scalar backwardH = forwardH;
    auto forwardJumpPosition = jumps.begin();
    auto backwardJumpPosition = jumps.rbegin();
    auto nextForwardJump = [&]() {
        return forwardJumpPosition == jumps.end() ? max : *forwardJumpPosition;
    };
    auto nextBackwardJump = [&]() {
        return backwardJumpPosition == jumps.rend() ? min : *backwardJumpPosition;
    };
    forward.emplace_back(
            automaticNextSector(problem, forwardH, min, std::min(mid, nextForwardJump()), tolerance,
                                Direction::forward));
    backward.emplace_back(
            automaticNextSector(problem, backwardH, std::max(mid, nextBackwardJump()), max, tolerance,
                                Direction::backward));


    while (forward.back()->max != backward.back()->min) {
        if (Sector::compare(*backward.back(), *forward.back())) {
            if (forward.back()->max == nextForwardJump())
                ++forwardJumpPosition;
            forward.emplace_back(automaticNextSector(
                    problem, forwardH,
                    forward.back()->max, std::min(backward.back()->min, nextForwardJump()),
                    tolerance, Direction::forward));
        } else {
            if (backward.back()->min == nextBackwardJump())
                ++backwardJumpPosition;
            backward.emplace_back(automaticNextSector(
                    problem, backwardH,
                    std::max(forward.back()->max, nextBackwardJump()), backward.back()->min,
                    tolerance, Direction::backward));
        }
    }

    int matchIndex = forward.size() - 1;
    forward.insert(forward.end(), make_move_iterator(backward.rbegin()), make_move_iterator(backward.rend()));
    return {std::move(forward), matchIndex};
};

#define INSTANTIATE_SECTOR_BUILDER(Problem) \
namespace matslise::sector_builder {                  \
    template class AutomaticSectorBuilder<Problem>; \
    template class UniformSectorBuilder<Problem>; \
}

#define INSTANTIATE_MATSLISE(Scalar) \
INSTANTIATE_SECTOR_BUILDER(Matslise<Scalar>)

#ifdef WITH_MATSCS

#include "../matscs.h"

#define INSTANTIATE_MORE(Scalar) \
INSTANTIATE_SECTOR_BUILDER(Matscs<Scalar>)
#endif


#include "instantiate.h"
