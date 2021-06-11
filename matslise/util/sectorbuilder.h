#ifndef MATSLISE_SECTORBUILDER_H
#define MATSLISE_SECTORBUILDER_H

#include <vector>
#include <optional>
#include <type_traits>
#include <functional>
#include "value_ptr.h"

namespace matslise {
    template<typename Problem>
    struct SectorBuilderReturn {
        std::vector<matslise::value_ptr<typename Problem::Sector>> sectors;
        // The index of the sector left of the matching point
        int matchIndex;
    };

    template<typename Problem, typename Scalar=typename Problem::Scalar>
    using SectorBuilder = std::function<SectorBuilderReturn<Problem>(
            const Problem *, const Scalar &min, const Scalar &max)>;

    namespace sector_builder {
        template<typename Problem>
        SectorBuilder<Problem> uniform(int sectorCount);

        template<typename Problem, bool parallel = false>
        SectorBuilder<Problem> automatic(const typename Problem::Scalar &tolerance);

        template<typename Problem, bool parallel = false>
        SectorBuilder<Problem> getOrAutomatic(const std::optional<SectorBuilder<Problem>> &sectorBuilder,
                                              const typename Problem::Scalar &tolerance) {
            if (sectorBuilder.has_value())
                return sectorBuilder.value();
            return automatic<Problem, parallel>(tolerance);
        }
    }
}

#endif //MATSLISE_SECTORBUILDER_H
