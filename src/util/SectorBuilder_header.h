#ifndef MATSLISE_SECTORBUILDER_HEADER_H
#define MATSLISE_SECTORBUILDER_HEADER_H

#include <functional>

namespace matslise {
    template<typename Problem>
    class SectorBuilder {
    public:
        virtual void build(Problem *, typename Problem::Scalar min, typename Problem::Scalar max) const = 0;

        virtual ~SectorBuilder() = default;;
    };

    namespace sectorbuilder {
        template<typename Problem>
        inline bool compareSectors(const typename Problem::Sector &a, const typename Problem::Sector &b);

        template<typename Problem>
        class Uniform : public SectorBuilder<Problem> {
        public:
            int sectorCount;

            explicit Uniform(int sectorCount) : sectorCount(sectorCount) {

            }

            virtual void build(Problem *ms, typename Problem::Scalar min, typename Problem::Scalar max) const;
        };

        template<typename Problem>
        class Auto : public SectorBuilder<Problem> {
        public:
            typename Problem::Scalar tol;

            explicit Auto(typename Problem::Scalar tol) : tol(tol) {

            }

            void build(Problem *ms, typename Problem::Scalar min, typename Problem::Scalar max) const override;

            template<bool forward>
            typename Problem::Sector *nextSector(
                    const Problem *ms, const typename Problem::Scalar &h, const typename Problem::Scalar &left,
                    const typename Problem::Scalar &right) const;
        };

        template<typename Problem>
        class Custom : public SectorBuilder<Problem> {
            std::function<void(Problem *, typename Problem::Scalar, typename Problem::Scalar)> builder;
        public:
            explicit Custom(std::function<void(Problem *, typename Problem::Scalar, typename Problem::Scalar)> builder)
                    : builder(builder) {}

            virtual void build(Problem *p, typename Problem::Scalar min, typename Problem::Scalar max) const {
                builder(p, min, max);
            }
        };
    }
}

#endif //MATSLISE_SECTORBUILDER_HEADER_H
