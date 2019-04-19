//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_MATSLISE_H
#define SCHRODINGER_MATSLISE_H

#include <cmath>
#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <functional>
#include "Array2D.h"
#include "util/eigen.h"
#include "util/y.h"
#include "Evaluator.h"

#define MATSLISE_HMAX_delta 15
#define MATSLISE_ETA_delta 7
#define MATSLISE_ETA_h 9
#define MATSLISE_N 15

using namespace Eigen;

namespace matslise {
    class Matslise;

    namespace matslise_util {
        class Sector;

        class SectorBuilderUniform {
        public:
            int sectorCount;

            SectorBuilderUniform(int sectorCount) : sectorCount(sectorCount) {

            }

            void build(Matslise *ms) const;
        };

        class SectorBuilderAuto {
        public:
            double tol;

            SectorBuilderAuto(double tol) : tol(tol) {

            }

            void build(Matslise *ms) const;

        private:
            template<bool forward>
            Sector *nextSector(Matslise *ms, double h, double left, double right) const;
        };

        class SectorBuilder {
        private:
            union {
                SectorBuilderAuto aut;
                SectorBuilderUniform uni;
            };
            enum {
                AUTO, UNIFORM
            } type;
        public:
            SectorBuilder(SectorBuilderAuto aut) {
                this->aut = aut;
                type = AUTO;
            }

            SectorBuilder(SectorBuilderUniform uni) {
                this->uni = uni;
                type = UNIFORM;
            }

            void build(Matslise *ms) const {
                switch (type) {
                    case UNIFORM:
                        return uni.build(ms);
                    case AUTO:
                        return aut.build(ms);
                }
            }
        };
    }

    class Matslise {
    public:
        std::function<double(double)> V;
        double xmin, xmax;
        int sectorCount;
        double match;
        matslise::matslise_util::Sector **sectors;
    public:
        static matslise_util::SectorBuilder UNIFORM(int sectorCount) {
            return matslise_util::SectorBuilder(matslise_util::SectorBuilderUniform(sectorCount));
        }

        static matslise_util::SectorBuilder AUTO(double tolerance) {
            return matslise_util::SectorBuilder(matslise_util::SectorBuilderAuto(tolerance));
        }

        Matslise(std::function<double(double)> V, double xmin, double xmax,
                 const matslise_util::SectorBuilder &sectorBuilder) : V(V), xmin(xmin), xmax(xmax) {
            sectorBuilder.build(this);
        }

        Matslise(std::function<double(double)> V, double xmin, double xmax, int sectorCount)
                : Matslise(V, xmin, xmax, UNIFORM(sectorCount)) {};

        std::pair<matslise::Y<>, double>
        propagate(double E, const matslise::Y<> &y, double a, double b, bool use_h = true) const;

        Eigen::Array<matslise::Y<>, Eigen::Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &left, const matslise::Y<> &right,
                             const Eigen::ArrayXd &x) const;

        std::tuple<double, double, double>
        calculateError(double E, const matslise::Y<> &left, const matslise::Y<> &right, bool use_h = true) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &left,
                           const matslise::Y<> &right) const;

        double computeEigenvalueError(double E, const matslise::Y<> &left, const matslise::Y<> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &left,
                                  const matslise::Y<> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, int Imin, int Imax,
                           const matslise::Y<> &left,
                           const matslise::Y<> &right) const;

        std::function<Y<>(double)> eigenfunctionCalculator(
                double E, const matslise::Y<> &left, const matslise::Y<> &right) const;

        virtual ~Matslise();
    };

    class HalfRange {
    public:
        const Matslise *ms;
        static const int AUTO = -1;
        static const int ODD = 0;
        static const int EVEN = 1;
    public:
        HalfRange(std::function<double(double)> V, double xmax,
                  const matslise::matslise_util::SectorBuilder &sectorBuilder);

        Array<matslise::Y<>, Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &side, const ArrayXd &x, int even = AUTO) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &side) const;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &side) const;

        double computeEigenvalueError(double E, const matslise::Y<> &side, int even = AUTO) const;

        std::function<Y<>(double)> eigenfunctionCalculator(double E, const matslise::Y<> &left, int even = AUTO) const;

        virtual ~HalfRange();
    };

    namespace matslise_util {
        class Sector {
        public:
            Matslise *s;
            double *vs;
            double xmin, xmax, h;
            Array2D<Matrix2d, MATSLISE_ETA_delta, MATSLISE_HMAX_delta>
                    t_coeff;
            Matrix2d t_coeff_h[MATSLISE_ETA_h];

            Sector(Matslise *problem, double xmin, double xmax);

            void calculateTCoeffs();

            T<> calculateT(double E, bool use_h = true) const;

            T<> calculateT(double E, double delta, bool use_h = true) const;

            Y<> propagate(double E, const Y<> &y0, bool forward, bool use_h = true) const;

            Y<> propagate(double E, const Y<> &y0, double delta, bool use_h = true) const;

            Y<> propagate(double E, const Y<> &y0, double delta, double &theta, bool use_h = true) const;

            double prufer(double E, double delta, const Y<> &y0, const Y<> &y1) const;

            double calculateError() const;

            virtual ~Sector();
        };
    };
}


#endif //SCHRODINGER_MATSLISE_H