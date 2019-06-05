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
#include "util/SectorBuilder.h"
#include <memory>

#define MATSLISE_HMAX_delta 15
#define MATSLISE_ETA_delta 7
#define MATSLISE_ETA_h 9
#define MATSLISE_N 15

using namespace Eigen;

namespace matslise {
    class AbstractMatslise {
    public:
        virtual std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &left,
                           const matslise::Y<> &right) const = 0;

        virtual std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &left) const {
            return computeEigenvalues(Emin, Emax, left, left);
        }

        virtual double
        computeEigenvalueError(double E, const matslise::Y<> &left, const matslise::Y<> &right) const = 0;

        virtual double computeEigenvalueError(double E, const matslise::Y<> &left) const {
            return computeEigenvalueError(E, left, left);
        };

        virtual std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &left,
                                  const matslise::Y<> &right) const = 0;

        virtual std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &left) const {
            return computeEigenvaluesByIndex(Imin, Imax, left, left);
        };


        virtual Eigen::Array<matslise::Y<>, Eigen::Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &left, const matslise::Y<> &right,
                             const Eigen::ArrayXd &x) const = 0;

        virtual Eigen::Array<matslise::Y<>, Eigen::Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &left, const Eigen::ArrayXd &x) const {
            return computeEigenfunction(E, left, left, x);
        }

        virtual std::function<Y<>(double)> eigenfunctionCalculator(
                double E, const matslise::Y<> &left, const matslise::Y<> &right) const = 0;

        virtual std::function<Y<>(double)> eigenfunctionCalculator(double E, const matslise::Y<> &left) const {
            return eigenfunctionCalculator(E, left, left);
        };

        virtual ~AbstractMatslise() {}
    };

    class Matslise : public AbstractMatslise {
    public:
        class Sector;

        std::function<double(double)> V;
        double xmin, xmax;
        int sectorCount;
        double match;
        Matslise::Sector **sectors;
    public:
        static std::shared_ptr<matslise::SectorBuilder<Matslise>> UNIFORM(int sectorCount) {
            return std::shared_ptr<matslise::SectorBuilder<Matslise>>(
                    new matslise::sectorbuilder::Uniform<Matslise>(sectorCount));
        }

        static std::shared_ptr<matslise::SectorBuilder<Matslise>> AUTO(double tolerance) {
            return std::shared_ptr<matslise::SectorBuilder<Matslise>>(
                    new matslise::sectorbuilder::Auto<Matslise>(tolerance));
        }

        Matslise(std::function<double(double)> V, double xmin, double xmax,
                 std::shared_ptr<matslise::SectorBuilder<Matslise>> sectorBuilder) : V(V), xmin(xmin), xmax(xmax) {
            sectorBuilder->build(this, xmin, xmax);
        }

        Matslise(std::function<double(double)> V, double xmin, double xmax, int sectorCount)
                : Matslise(V, xmin, xmax, UNIFORM(sectorCount)) {};

        bool contains(double point) const {
            return point <= xmax && point >= xmin;
        }

        std::pair<matslise::Y<>, double>
        propagate(double E, const matslise::Y<> &y, double a, double b, bool use_h = true) const;

        std::tuple<double, double, double>
        calculateError(double E, const matslise::Y<> &left, const matslise::Y<> &right, bool use_h = true) const;

    public: // Override
        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &left,
                           const matslise::Y<> &right) const override;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &left,
                                  const matslise::Y<> &right) const override;

        double computeEigenvalueError(double E, const matslise::Y<> &left, const matslise::Y<> &right) const override;

        Eigen::Array<matslise::Y<>, Eigen::Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &left, const matslise::Y<> &right,
                             const Eigen::ArrayXd &x) const override;

        std::function<Y<>(double)>
        eigenfunctionCalculator(double E, const matslise::Y<> &left, const matslise::Y<> &right) const override;

        virtual ~Matslise();

    public:
        class Sector {
        public:
            Array2D<Matrix2d, MATSLISE_ETA_delta, MATSLISE_HMAX_delta>
                    t_coeff;
            Matrix2d t_coeff_h[MATSLISE_ETA_h];
            Matslise *s;
            double *vs;
            double min, max, h;
            bool backward;

            Sector(Matslise *problem, double min, double max, bool backward);

            void calculateTCoeffs();

            bool contains(double point) const {
                return point <= max && point >= min;
            }

            T<> calculateT(double E, bool use_h = true) const;

            T<> calculateT(double E, double delta, bool use_h = true) const;

            matslise::Y<> propagate(double E, const matslise::Y<> &y, bool forward, bool use_h = true) const;

            Y<> propagate(double E, const Y<> &y0, double a, double b, double &theta, bool use_h = true) const;

            double prufer(double E, double delta, const Y<> &y0, const Y<> &y1) const;

            double calculateError() const;

            virtual ~Sector();
        };
    };

    class HalfRange : public AbstractMatslise {
    public:
        static void checkSymmetry(const matslise::Y<> &left, const matslise::Y<> &right) {
            if (left != right)
                throw std::runtime_error("Halfrange::checkSymmetry(), left and right sides have to be identical");
        }

    public:
        const Matslise *ms;
        static const int AUTO = -1;
        static const int ODD = 0;
        static const int EVEN = 1;
    public:
        HalfRange(std::function<double(double)> V, double xmax,
                  std::shared_ptr<matslise::SectorBuilder<Matslise>> sectorBuilder);

        using AbstractMatslise::computeEigenvalues;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<> &left,
                           const matslise::Y<> &right) const override;

        using AbstractMatslise::computeEigenvaluesByIndex;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<> &left,
                                  const matslise::Y<> &right) const override;


        using AbstractMatslise::computeEigenvalueError;

        double
        computeEigenvalueError(double E, const matslise::Y<> &side, int even) const;

        double computeEigenvalueError(double E, const matslise::Y<> &left, const matslise::Y<> &right) const override {
            checkSymmetry(left, right);
            return computeEigenvalueError(E, left, AUTO);
        }

        using AbstractMatslise::computeEigenfunction;

        Array<matslise::Y<>, Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &left, const matslise::Y<> &right,
                             const ArrayXd &x) const override {
            checkSymmetry(left, right);
            return computeEigenfunction(E, left, x, AUTO);
        }

        Array<matslise::Y<>, Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<> &side, const ArrayXd &x, int even) const;

        using AbstractMatslise::eigenfunctionCalculator;

        std::function<Y<>(double)>
        eigenfunctionCalculator(double E, const matslise::Y<> &left, const matslise::Y<> &right) const override {
            checkSymmetry(left, right);
            return eigenfunctionCalculator(E, left, AUTO);
        }

        std::function<Y<>(double)>
        eigenfunctionCalculator(double E, const matslise::Y<> &side, int even) const;

        virtual ~HalfRange();
    };
}


#endif //SCHRODINGER_MATSLISE_H