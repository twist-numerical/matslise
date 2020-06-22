#ifndef MATSLISE_SECTORBUILDER_H
#define MATSLISE_SECTORBUILDER_H

#include <vector>
#include <type_traits>
#include <functional>

namespace matslise {
    template<typename Problem>
    struct SectorBuilderReturn {
        std::vector<typename Problem::Sector *> sectors;
        // The index of the sector left of the matching point
        int matchIndex;
    };

    template<typename Problem, typename Scalar=typename Problem::Scalar>
    using SectorBuilder = std::function<SectorBuilderReturn<Problem>(
            const Problem *, const Scalar &min, const Scalar &max)>;

    namespace sector_builder {
        template<typename Problem>
        SectorBuilder<Problem> uniform(int sectorCount);

        template<typename Problem>
        SectorBuilder<Problem> automatic(const typename Problem::Scalar &tolerance);
    }
}

#endif //MATSLISE_SECTORBUILDER_H
