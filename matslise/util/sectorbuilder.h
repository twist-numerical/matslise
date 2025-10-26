#ifndef MATSLISE_SECTORBUILDER_H
#define MATSLISE_SECTORBUILDER_H

#include <vector>
#include <memory>

namespace matslise::sector_builder {
    template<typename Problem>
    struct SectorBuilderReturn {
        std::vector<std::unique_ptr<typename Problem::Sector>> sectors{};
        // The index of the sector left of the matching point
        int matchIndex{};
    };

    template<typename Problem, typename Scalar=typename Problem::Scalar>
    class SectorBuilder {
    public:
        virtual ~SectorBuilder() = default;
        virtual SectorBuilderReturn<Problem>
        operator()(const Problem *, const Scalar &min, const Scalar &max) const = 0;
    };

    template<typename Problem, typename Scalar=typename Problem::Scalar>
    class AutomaticSectorBuilder : public SectorBuilder<Problem, Scalar> {
    public:
        std::vector<Scalar> jumps{};
        Scalar tolerance;

        AutomaticSectorBuilder(const Scalar &tolerance_) : tolerance(tolerance_) {
        }

        SectorBuilderReturn<Problem> operator()(const Problem *, const Scalar &min, const Scalar &max) const override;
    };

    template<typename Problem, typename Scalar=typename Problem::Scalar>
    class UniformSectorBuilder : public SectorBuilder<Problem, Scalar> {
    public:
        int sectorCount;

        UniformSectorBuilder(int sectorCount_) : sectorCount(sectorCount_) {
        }

        SectorBuilderReturn<Problem> operator()(const Problem *, const Scalar &min, const Scalar &max) const override;
    };
}

#endif //MATSLISE_SECTORBUILDER_H
