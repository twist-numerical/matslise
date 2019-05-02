//
// Created by toon on 4/23/19.
//

#ifndef MATSLISE_SECTORBUILDER_H
#define MATSLISE_SECTORBUILDER_H

#include <functional>

namespace matslise {
    template<typename Problem>
    class SectorBuilder {
    public:
        virtual void build(Problem *, double min, double max) const = 0;

        virtual ~SectorBuilder() {};
    };

    namespace sectorbuilder {

        template<typename Problem>
        class Uniform : public SectorBuilder<Problem> {
        public:
            int sectorCount;

            Uniform(int sectorCount) : sectorCount(sectorCount) {

            }

            virtual void build(Problem *, double min, double max) const;
        };

        template<typename Problem>
        class Auto : public SectorBuilder<Problem> {
        public:
            double tol;

            Auto(double tol) : tol(tol) {

            }

            virtual void build(Problem *, double min, double max) const;

        private:
            template<bool forward>
            typename Problem::Sector *nextSector(Problem *ms, double h, double left, double right) const;
        };

        template<typename Problem>
        class Custom : public SectorBuilder<Problem> {
            std::function<void(Problem *, double, double)> builder;
        public:
            Custom(std::function<void(Problem *, double, double)> builder) : builder(builder) {}

            virtual void build(Problem *p, double min, double max) const {
                builder(p, min, max);
            }
        };
    }
}

#endif //MATSLISE_SECTORBUILDER_H
