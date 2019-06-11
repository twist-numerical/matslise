//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_MATSCS_H
#define SCHRODINGER_MATSCS_H

#include <ostream>
#include <array>
#include <vector>
#include <functional>

#define MATSCS_HMAX_delta 7
#define MATSCS_ETA_delta 3
#define MATSCS_ETA_h 6
#define MATSCS_N 9

namespace matslise {
    template<typename _Scalar>
    class Matscs {
    public:
        typedef _Scalar Scalar;

        class Sector;

        std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(Scalar)> V;
        int n;
        Scalar xmin, xmax;
        int sectorCount;
        Matscs::Sector **sectors;
        Scalar match;
    public:
        static std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>> UNIFORM(int sectorCount) {
            return std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>>(
                    new matslise::sectorbuilder::Uniform<Matscs<Scalar>>(sectorCount));
        }

        static std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>> AUTO(Scalar tolerance) {
            return std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>>(
                    new matslise::sectorbuilder::Auto<Matscs<Scalar>>(tolerance));
        }

        Matscs(std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(Scalar)> V, int n, Scalar xmin,
               Scalar xmax,
               std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>> sectorBuilder) : V(V), n(n), xmin(xmin),
                                                                                         xmax(xmax) {
            sectorBuilder->build(this, xmin, xmax);
        }

        Matscs(std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(Scalar)> V, int n, Scalar xmin,
               Scalar xmax,
               int sectorCount)
                : Matscs(V, n, xmin, xmax, UNIFORM(sectorCount)) {};

        bool contains(Scalar point) const {
            return point <= xmax && point >= xmin;
        }

        template<int r>
        matslise::Y<Scalar, Eigen::Dynamic, r>
        propagate(Scalar E, const matslise::Y<Scalar, Eigen::Dynamic, r> &y, Scalar a, Scalar b,
                  bool use_h = true) const;

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
        propagatePsi(Scalar E, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &psi, Scalar a,
                     Scalar b) const;

        std::vector<matslise::Y<Scalar, Eigen::Dynamic>> *computeEigenfunction(Scalar E, std::vector<Scalar> &x) const;

        std::function<Y<Scalar, Eigen::Dynamic, 1>(Scalar)> eigenfunctionCalculator(
                Scalar E, const matslise::Y<Scalar, Eigen::Dynamic, 1> &left,
                const matslise::Y<Scalar, Eigen::Dynamic, 1> &right);

        ~Matscs();

        class Sector {
        private:
            const Matscs *s;
        public:
            Eigen::Array<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, MATSCS_ETA_delta, MATSCS_HMAX_delta>
                    t_coeff;
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> t_coeff_h[MATSCS_ETA_h];
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> diagonalize;
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vs;
            Scalar min, max, h;
            bool backward;

            Sector(const Matscs *problem, Scalar min, Scalar max, bool backward = false);

            void calculateTCoeffs();

            bool contains(Scalar point) const {
                return point <= max && point >= min;
            }

            T<Scalar, Eigen::Dynamic> calculateT(Scalar E, bool use_h = true) const;

            T<Scalar, Eigen::Dynamic> calculateT(Scalar E, Scalar delta, bool use_h = true) const;

            template<int r = Eigen::Dynamic>
            Y<Scalar, Eigen::Dynamic, r>
            propagate(Scalar E, const Y<Scalar, Eigen::Dynamic, r> &y0, Scalar a, Scalar b, bool use_h = true) const;

            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
            propagatePsi(Scalar E, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &psi, Scalar a,
                         Scalar b) const;

            Scalar calculateError() const;

            ~Sector();

        };
    };
}

#include "util/SectorBuilder.h"

#endif