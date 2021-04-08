#include "./sectorbuilder.h"
#include "../matslise.h"
#include "./heap.h"
#include <algorithm>
#include <vector>
#include <tuple>


using namespace matslise;

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

template<typename Problem>
SectorBuilder<Problem> sector_builder::uniform(int sectorCount) {
    typedef typename Problem::Scalar Scalar;
    typedef typename Problem::Sector Sector;
    return [sectorCount](const Problem *problem, const Scalar &min, const Scalar &max)
            -> SectorBuilderReturn<Problem> {
        Scalar h = (max - min) / sectorCount;

        if (sectorCount == 1) {
            return {{value_ptr<Sector>(new Sector(problem, min, max, forward))}, 0};
        }

        std::vector<value_ptr<Sector>> sectors;
        sectors.resize(sectorCount);
#ifdef MATSLISE_parallel
#pragma omp parallel for shared(min, h)
        for (int i = 0; i < sectorCount; ++i) {
            sectors[i].reset(new Sector(problem, min + i * h, max - (sectorCount - i - 1) * h, none));
        }
        int matchIndex = 0;
        for (int i = 1; i < sectorCount; ++i) {
            if (Sector::compare(*sectors[i], *sectors[matchIndex]))
                matchIndex = i - 1;
        }

#pragma omp parallel for
        for (int i = 0; i < sectorCount; ++i) {
            sectors[i]->setDirection(i <= matchIndex ? forward : backward);
        }

#else
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
        h *= std::max(Scalar(.1), .99 * Scalar(pow(tolerance / error, Scalar(1. / (Problem::order - 1)))));
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
SectorBuilder<Problem> automaticSequential(const typename Problem::Scalar &tolerance) {
    typedef typename Problem::Scalar Scalar;
    typedef typename Problem::Sector Sector;
    return [tolerance](const Problem *problem, const Scalar &min, const Scalar &max)
            -> SectorBuilderReturn<Problem> {
        std::vector<value_ptr<typename Problem::Sector>> forward;
        std::vector<value_ptr<typename Problem::Sector>> backward;
        // It shouldn't be the true middle
        Scalar mid = 0.4956864123 * max + 0.5043135877 * min;
        Scalar forwardH = .33 * (max - min);
        Scalar backwardH = forwardH;
        forward.emplace_back(automaticNextSector(problem, forwardH, min, mid, tolerance, Direction::forward));
        backward.emplace_back(automaticNextSector(problem, backwardH, mid, max, tolerance, Direction::backward));


        while (forward.back()->max != backward.back()->min) {
            if (Sector::compare(*backward.back(), *forward.back())) {
                forward.emplace_back(automaticNextSector(
                        problem, forwardH, forward.back()->max,
                        backward.back()->min, tolerance, Direction::forward));
            } else {
                backward.emplace_back(automaticNextSector(
                        problem, backwardH, forward.back()->max,
                        backward.back()->min, tolerance, Direction::backward));
            }
        }

        int matchIndex = forward.size() - 1;
        forward.insert(forward.end(), make_move_iterator(backward.rbegin()), make_move_iterator(backward.rend()));
        return {std::move(forward), matchIndex};
    };
};

#ifdef MATSLISE_parallel

template<typename Problem>
SectorBuilder<Problem> automaticParallel(const typename Problem::Scalar &tolerance) {
    typedef typename Problem::Scalar Scalar;
    typedef typename Problem::Sector Sector;
    typedef std::tuple<Scalar, Sector *, int> Item;

    return [tolerance](const Problem *problem, const Scalar &min, const Scalar &max)
            -> SectorBuilderReturn<Problem> {
        Heap<Item> heap([](const Item &a, const Item &b) { return std::get<0>(a) < std::get<0>(b); });
        Scalar maxError = 100 * tolerance;
        heap.emplace(maxError, new Sector(problem, min, max, none), 0);

        int maxDepth = 20;
#pragma omp parallel
        {
#pragma omp single nowait
            {
                while (maxError > tolerance) {
                    bool empty = false;
                    while (!empty && maxError > tolerance) {
                        Sector *sector;
                        int depth;
#pragma omp critical (edit_heap)
                        {
                            sector = std::get<1>(heap.front());
                            depth = std::get<2>(heap.front());
                            heap.pop();
                        }
                        Scalar m1 = 0.499 * sector->max + 0.501 * sector->min;
                        // Scalar m1 = 0.325 * sector->max + 0.675 * sector->min;
                        // Scalar m2 = 0.675 * sector->max + 0.325 * sector->min;
                        for (auto &bounds : (std::pair<Scalar, Scalar>[]) {
                                {sector->min, m1},
                                {m1,          sector->max}
                                // {m1,          m2},
                                // {m2,          sector->max}
                        }) {
#pragma omp task
                            {

                                Sector *newSector = refineSector(*sector, problem, bounds.first, bounds.second,
                                                                 forward);
                                Scalar newError = depth < maxDepth ? newSector->error() : 0;
                                if (!(newError < 100 * tolerance)) { // (! . < .) also correct for nan
                                    newError = 100 * tolerance;
                                }
#pragma omp critical (edit_heap)
                                {
                                    heap.emplace(newError, newSector, depth + 1);
                                }
                            }
                        }
#pragma omp critical (edit_heap)
                        {
                            empty = heap.empty();
                            maxError = std::get<0>(heap.front());
                        }
                    }
#pragma omp taskwait
                }
            }
        }

        SectorBuilderReturn<Problem> result;
        result.sectors.reserve(heap.size());
        for (auto &item : heap.data())
            result.sectors.template emplace_back().reset(std::get<1>(item));
        std::sort(result.sectors.begin(), result.sectors.end(), [](const value_ptr<Sector> &a, const value_ptr<Sector> &b) {
            return a->min < b->min;
        });

        auto match = result.sectors.begin();
        for (auto i = match + 1; i != result.sectors.end(); ++i) {
            if (Sector::compare(**i, **match))
                match = i - 1;
        }
        result.matchIndex = std::distance(result.sectors.begin(), match);

#pragma omp parallel for
        for (int i = 0; i < result.matchIndex; ++i) {
            result.sectors[i]->setDirection(i <= result.matchIndex ? forward : backward);
        }

        return result;
    };
}

#endif

template<typename Problem, bool parallel>
SectorBuilder<Problem> sector_builder::automatic(const typename Problem::Scalar &tolerance) {
#ifdef MATSLISE_parallel
    if constexpr(parallel) {
        return automaticParallel<Problem>(tolerance);
    } else {
        return automaticSequential<Problem>(tolerance);
    }
#else
    return automaticSequential<Problem>(tolerance);
#endif
}

#define INSTANTIATE_SECTOR_BUILDER(Problem) \
template SectorBuilder<Problem> sector_builder::uniform<Problem>(int); \
template SectorBuilder<Problem> sector_builder::automatic<Problem, false>(const Problem::Scalar&); \
template SectorBuilder<Problem> sector_builder::automatic<Problem, true>(const Problem::Scalar&);

#define INSTANTIATE_MORE(Scalar) \
INSTANTIATE_SECTOR_BUILDER(Matslise<Scalar>) \
INSTANTIATE_SECTOR_BUILDER(Matscs<Scalar>) \
INSTANTIATE_SECTOR_BUILDER(Matslise2D<Scalar>) \
INSTANTIATE_SECTOR_BUILDER(Matslise3D<Scalar>)

#include "./instantiate.h"
