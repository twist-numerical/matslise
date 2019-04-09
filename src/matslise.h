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

    namespace matslise_util {
        class Sector;

        class EigenfunctionCalculator;
    }

    class Matslise {
    public:
        std::function<double(double)> V;
        double xmin, xmax;
        int sectorCount;
        double match;
        matslise::matslise_util::Sector **sectors;
    public:
        Matslise(std::function<double(double)> V, double xmin, double xmax, int sectorCount);

        std::pair<matslise::Y<>, double>
        propagate(double E, const matslise::Y<> &y, double a, double b) const;

        Eigen::Array<matslise::Y<>, Eigen::Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &left, const matslise::Y<> &right,
                             const Eigen::ArrayXd &x) const;

        std::tuple<double, double, double>
        calculateError(double E, const matslise::Y<> &left, const matslise::Y<> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &left,
                           const matslise::Y<> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &left,
                                  const matslise::Y<> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, int Imin, int Imax,
                           const matslise::Y<> &left,
                           const matslise::Y<> &right) const;

        matslise::matslise_util::EigenfunctionCalculator *eigenfunctionCalculator(
                double E, const matslise::Y<> &left, const matslise::Y<> &right);

        virtual ~Matslise();
    };

    class HalfRange {
    public:
        const Matslise *ms;

    public:
        HalfRange(std::function<double(double)> V, double xmax, int sectorCount);

        Array<matslise::Y<>, Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &side, const ArrayXd &x) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &side) const;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &side) const;

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

            T<> calculateT(double E) const;

            T<> calculateT(double E, double delta) const;

            Y<> propagate(double E, const Y<> &y0, bool forward) const;

            Y<> propagate(double E, const Y<> &y0, double delta) const;

            Y<> propagate(double E, const Y<> &y0, double delta, double &theta) const;

            double prufer(double E, double delta, const Y<> &y0, const Y<> &y1) const;

            virtual ~Sector();
        };

        class EigenfunctionCalculator : public Evaluator<Y<>, double> {
        public:
            Matslise *ms;
            double E;
            std::vector<Y<>> ys;

            EigenfunctionCalculator(Matslise *ms, double E, const Y<> &left, const Y<> &right);

            virtual Y<> eval(double x) const;
        };
    };
}


#endif //SCHRODINGER_MATSLISE_H