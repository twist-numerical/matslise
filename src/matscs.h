//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_MATSCS_H
#define SCHRODINGER_MATSCS_H

#include "matslise.h"
#include <ostream>
#include <array>
#include <vector>
#include <functional>
#include "Array2D.h"
#include "util/eigen.h"

#define MATSCS_HMAX_delta 7
#define MATSCS_ETA_delta 3
#define MATSCS_ETA_h 6
#define MATSCS_N 9

using namespace Eigen;


namespace matslise {
    class Matscs {
    public:
        class Sector;

        std::function<MatrixXd(double)> V;
        int n;
        double xmin, xmax;
        int sectorCount;
        Matscs::Sector **sectors;
        double match;
    public:
        static std::unique_ptr<matslise::SectorBuilder<Matscs>> UNIFORM(int sectorCount) {
            return std::unique_ptr<matslise::SectorBuilder<Matscs>>(
                    new matslise::sectorbuilder::Uniform<Matscs>(sectorCount));
        }

        static std::unique_ptr<matslise::SectorBuilder<Matscs>> AUTO(double tolerance) {
            return std::unique_ptr<matslise::SectorBuilder<Matscs>>(
                    new matslise::sectorbuilder::Auto<Matscs>(tolerance));
        }

        Matscs(std::function<MatrixXd(double)> V, int n, double xmin, double xmax,
               std::unique_ptr<matslise::SectorBuilder<Matscs>> sectorBuilder) : V(V), n(n), xmin(xmin), xmax(xmax) {
            sectorBuilder->build(this, xmin, xmax);
        }

        Matscs(std::function<MatrixXd(double)> V, int n, double xmin, double xmax, int sectorCount)
                : Matscs(V, n, xmin, xmax, UNIFORM(sectorCount)) {};

        template<int r>
        matslise::Y<Eigen::Dynamic, r>
        propagate(double E, const matslise::Y<Eigen::Dynamic, r> &y, double a, double b, bool use_h = true) const;

        MatrixXd propagatePsi(double E, const MatrixXd &psi, double a, double b) const;

        std::vector<matslise::Y<Eigen::Dynamic>> *computeEigenfunction(double E, std::vector<double> &x) const;

        std::function<Y<Dynamic, 1>(double)> eigenfunctionCalculator(
                double E, const matslise::Y<Eigen::Dynamic, 1> &left, const matslise::Y<Eigen::Dynamic, 1> &right);

        ~Matscs();

        class Sector {
        private:
            const Matscs *s;
        public:
            MatrixXd *vs;
            double min, max, h;
            Array2D<MatrixXd, MATSCS_ETA_delta, MATSCS_HMAX_delta> t_coeff;
            MatrixXd t_coeff_h[MATSCS_ETA_h];
            MatrixXd D;

            Sector(const Matscs *problem, double min, double max);

            void calculateTCoeffs();

            T<Eigen::Dynamic> calculateT(double E, bool use_h = true) const;

            T<Eigen::Dynamic> calculateT(double E, double delta, bool use_h = true) const;

            template<int r = Eigen::Dynamic>
            Y<Eigen::Dynamic, r>
            propagate(double E, const Y<Eigen::Dynamic, r> &y0, double delta, bool use_h = true) const;

            MatrixXd propagatePsi(double E, const MatrixXd &psi, double delta) const;

            double calculateError() const;

            ~Sector();

        };
    };
}

#endif