//
// Created by toon on 6/13/18.
//

#ifndef MATSLISE_SE2D_H
#define MATSLISE_SE2D_H


#include <functional>
#include <vector>
#include "matslise.h"
#include "matscs.h"

namespace matslise {
    namespace SEnD_util {
        const std::function<bool(std::pair<double, double>, std::pair<double, double>)>
                NEWTON_RAPHSON_SORTER =
                [](const std::pair<double, double> &a, const std::pair<double, double> &b) {
                    if (abs(a.first) > 100 || abs(b.first) > 100)
                        return abs(a.first) < abs(b.first);
                    return abs(a.first / a.second) < abs(b.first / b.second);
                };

        const std::function<bool(std::pair<double, double>, std::pair<double, double>)>
                ABS_SORTER =
                [](const std::pair<double, double> &a, const std::pair<double, double> &b) {
                    return abs(a.first) < abs(b.first);
                };

        template<int n>
        class Sector;

        template<int n>
        class EigenfunctionCalculator;
    }

    template<int n>
    struct Options {
        Options<n - 1> nestedOptions;
        int _sectorCount = 22;
        int _stepsPerSector = 3;
        int _N = 12;
        int _gridPoints = 50;

        Options<n> &nested(Options<n - 1> &nested) {
            nestedOptions = nested;
            return *this;
        }

        Options<n> &sectorCount(int count) {
            _sectorCount = count;
            return *this;
        }

        Options<n> &stepsPerSector(int count) {
            _stepsPerSector = count;
            return *this;
        }

        Options<n> &N(int value) {
            _N = value;
            return *this;
        }

        Options<n> &gridPoints(int value) {
            _gridPoints = value;
            return *this;
        }
    };

    template<>
    struct Options<1> {
        int _sectorCount = 26;

        Options<1> &sectorCount(int count) {
            _sectorCount = count;
            return *this;
        }
    };

    template<int n>
    struct dim {
    };

    template<>
    struct dim<2> {
        typedef std::function<double(double, double)> function;
    };

    template<>
    struct dim<3> {
        typedef std::function<double(double, double, double)> function;
    };

    template<int n>
    class Rectangle {
    public:
        Rectangle<n - 1> sub;
        double min, max;
    };

    template<>
    class Rectangle<0> {
    };

    template<int n = 2>
    class SEnD {
    public:
        typename dim<n>::function V;
        MatrixXd *M;
        Rectangle<n> domain;
        int sectorCount;
        SEnD_util::Sector<n> **sectors;
        ArrayXd xGrid;
        int N;

    public:
        SEnD(typename dim<n>::function V, const Rectangle<n> &domain, const Options<n> &options);

        std::pair<MatrixXd, MatrixXd> calculateErrorMatrix(double E) const;

        std::pair<double, double> calculateError(double E, const std::function<bool(
                std::pair<double, double>,
                std::pair<double, double>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER) const;

        std::vector<std::pair<double, double>> *calculateErrors(double E) const;

        std::vector<std::pair<double, double>> *sortedErrors(double E, const std::function<bool(
                std::pair<double, double>,
                std::pair<double, double>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER) const;

        std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

        std::vector<Eigen::ArrayXXd>
        computeEigenfunction(double E, const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const;

        Y<MatrixXd> *computeEigenfunctionSteps(double E) const;

        MatrixXd calculateM(int k) const;

        const SEnD_util::Sector<n> &getSector(double y) const;

        virtual ~SEnD();
    };

    namespace SEnD_util {
        template<int n = 2>
        class Sector {
        public:
            SEnD<n> *se2d;
            Matslise *matslise;
            Matscs *matscs;
            ArrayXd vbar;
            double ymin, ymax;

            double *eigenvalues;
            double *eigenfunctionsScaling;
            Eigen::ArrayXd *eigenfunctions;

            Sector(SEnD<n> *se2d, double ymin, double ymax, const Options<n> &options);

            matslise::Y<MatrixXd> propagate(double E, const matslise::Y<MatrixXd> &c, bool forward) const;

            matslise::Y<MatrixXd> propagate(double E, const matslise::Y<MatrixXd> &c, double y, bool forward) const;

            ArrayXd computeEigenfunction(int index, const ArrayXd &x) const;

            Eigen::MatrixXd calculateDeltaV(double y) const;

            virtual ~Sector();
        };


        template<int n = 2>
        class EigenfunctionCalculator : public Evaluator<double, double, double> {
        public:
            SEnD<n> *se2d;
            double E;
            std::vector<Y<VectorXd>> ys;

            EigenfunctionCalculator(SEnD<n> *se2d, double E);

            virtual double eval(double x, double y) const;
        };
    }
}

#endif //MATSLISE_SE2D_H
