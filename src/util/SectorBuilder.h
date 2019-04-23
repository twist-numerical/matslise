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
        virtual void build(Problem *) const = 0;
        virtual ~SectorBuilder() {};

    public:
        class Uniform : public SectorBuilder<Problem> {
        public:
            int sectorCount;

            Uniform(int sectorCount) : sectorCount(sectorCount) {

            }

            virtual void build(Problem *) const;
        };

        class Auto : public SectorBuilder<Problem> {
        public:
            double tol;

            Auto(double tol) : tol(tol) {

            }

            virtual void build(Problem *) const;

        private:
            template<bool forward>
            typename Problem::Sector *nextSector(Problem *ms, double h, double left, double right) const;
        };

        class Custom : public SectorBuilder<Problem> {
            std::function<typename Problem::Sector(Problem *)> builder;
        public:
            Custom(std::function<typename Problem::Sector(Problem *)> builder) : builder(builder) {}

            virtual void build(Problem *p) const {
                builder(p);
            }
        };
    };
}

#endif //MATSLISE_SECTORBUILDER_H
