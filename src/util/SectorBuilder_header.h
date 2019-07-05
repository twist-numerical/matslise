#ifndef MATSLISE_SECTORBUILDER_HEADER_H
#define MATSLISE_SECTORBUILDER_HEADER_H

#include <functional>

namespace matslise {
    template<typename Problem>
    class SectorBuilder {
    public:
        virtual void build(Problem *, typename Problem::Scalar min, typename Problem::Scalar max) const = 0;

        virtual ~SectorBuilder() {};
    };

    namespace sectorbuilder {
        template<typename Problem>
        inline bool compareSectors(typename Problem::Sector *a, typename Problem::Sector *b);

        template<typename Problem>
        class Uniform : public SectorBuilder<Problem> {
        public:
            int sectorCount;

            Uniform(int sectorCount) : sectorCount(sectorCount) {

            }

            virtual void build(Problem *ms, typename Problem::Scalar min, typename Problem::Scalar max) const;
        };

        template<typename Problem>
        class Auto : public SectorBuilder<Problem> {
        public:
            typename Problem::Scalar tol;

            Auto(typename Problem::Scalar tol) : tol(tol) {

            }

            virtual void build(Problem *ms, typename Problem::Scalar min, typename Problem::Scalar max) const override;

            template<bool forward>
            typename Problem::Sector *nextSector(Problem *ms, typename Problem::Scalar h, typename Problem::Scalar left,
                                                 typename Problem::Scalar right) const;
        };

        template<typename Problem>
        class Custom : public SectorBuilder<Problem> {
            std::function<void(Problem *, typename Problem::Scalar, typename Problem::Scalar)> builder;
        public:
            Custom(std::function<void(Problem *, typename Problem::Scalar, typename Problem::Scalar)> builder)
                    : builder(builder) {}

            virtual void build(Problem *p, typename Problem::Scalar min, typename Problem::Scalar max) const {
                builder(p, min, max);
            }
        };
    }
}

#endif //MATSLISE_SECTORBUILDER_HEADER_H
