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
#define MATSCS_ETA_h 4
#define MATSCS_N 5

using namespace Eigen;


namespace matslise {
    namespace matscs_util {
        class Sector;

        class EigenfunctionCalculator;
    }

    class Matscs {
    public:
        std::function<MatrixXd(double)> V;
        int n;
        double xmin, xmax;
        int sectorCount;
        matslise::matscs_util::Sector **sectors;
    public:
        Matscs(std::function<MatrixXd(double)> V, int n, double xmin, double xmax, int sectorCount);

        template<int r>
        matslise::Y<Eigen::Dynamic, r>
        propagate(double E, const matslise::Y<Eigen::Dynamic, r> &y, double a, double b) const;

        MatrixXd propagatePsi(double E, const MatrixXd &psi, double a, double b) const;

        std::vector<matslise::Y<Eigen::Dynamic>> *computeEigenfunction(double E, std::vector<double> &x) const;

        matslise::matscs_util::EigenfunctionCalculator *eigenfunctionCalculator(
                double E, const matslise::Y<Eigen::Dynamic, 1> &left, const matslise::Y<Eigen::Dynamic, 1> &right);

        virtual ~Matscs();
    };

    namespace matscs_util {
        class Sector {
        private:
            const Matscs *s;
        public:
            MatrixXd *vs;
            double xmin, xmax, h;
            Array2D<MatrixXd, MATSCS_ETA_delta, MATSCS_HMAX_delta> t_coeff;
            MatrixXd t_coeff_h[MATSCS_ETA_h];
            MatrixXd D;

            Sector(const Matscs *problem, double xmin, double xmax);

            void calculateTCoeffs();

            T<Eigen::Dynamic> calculateT(double E) const;

            T<Eigen::Dynamic> calculateT(double E, double delta) const;

            template<int r = Eigen::Dynamic>
            Y<Eigen::Dynamic, r> propagate(double E, const Y<Eigen::Dynamic, r> &y0, double delta) const;

            MatrixXd propagatePsi(double E, const MatrixXd &psi, double delta) const;

            virtual ~Sector();

        };


        class EigenfunctionCalculator : public Evaluator<Y<Eigen::Dynamic, 1>, double> {
        public:
            Matscs *ms;
            double E;
            std::vector<Y<Eigen::Dynamic, 1>> ys;

            EigenfunctionCalculator(Matscs *ms, double E, const Y<Eigen::Dynamic, 1> &left,
                                    const Y<Eigen::Dynamic, 1> &right);

            virtual Y<Eigen::Dynamic, 1> eval(double x) const;
        };
    }
}

#endif