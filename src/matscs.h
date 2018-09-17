//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_MATSCS_H
#define SCHRODINGER_MATSCS_H

#include <ostream>
#include <array>
#include <vector>
#include <functional>
#include "Array2D.h"
#include "matslise.h"
#include <Eigen/Dense>

#define MATSCS_HMAX 9
#define MATSCS_ETA 5

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
        MatrixXd zero, one;
    public:
        Matscs(std::function<MatrixXd(double)> V, int n, double xmin, double xmax, int sectorCount);

        template<typename Type, int... Args>
        matslise::Y<Matrix<Type, Args...>> propagate(double E, const matslise::Y<Matrix<Type, Args...>> &y, double a, double b) const;

        std::vector<matslise::Y<MatrixXd>> *computeEigenfunction(double E, std::vector<double> &x) const;

        matslise::matscs_util::EigenfunctionCalculator *eigenfunctionCalculator(
                double E, const matslise::Y<VectorXd> &left, const matslise::Y<VectorXd> &right);

        virtual ~Matscs();
    };

    namespace matscs_util {
        class Sector {
        private:
            const Matscs *s;
        public:
            MatrixXd *vs;
            double xmin, xmax, h;
            Array2D<Matrix2D<MatrixXd>, MATSCS_ETA, MATSCS_HMAX> t_coeff;
            Matrix2D<MatrixXd> t_coeff_h[MATSCS_ETA];
            MatrixXd D;

            Sector(const Matscs *problem, double xmin, double xmax);

            void calculateTCoeffs();

            T<MatrixXd> calculateT(double E) const;

            T<MatrixXd> calculateT(double E, double delta) const;

            template<typename Type, int... Args>
            Y<Matrix<Type, Args...>> propagate(double E, const Y<Matrix<Type, Args...>> &y0, double delta) const;

            virtual ~Sector();

        };


        class EigenfunctionCalculator : public Evaluator<Y<VectorXd>, double> {
        public:
            Matscs *ms;
            double E;
            std::vector<Y<VectorXd>> ys;

            EigenfunctionCalculator(Matscs *ms, double E, const Y<VectorXd> &left, const Y<VectorXd> &right);

            virtual Y<VectorXd> eval(double x) const;
        };
    }
}

#endif