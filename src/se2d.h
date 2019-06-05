//
// Created by toon on 6/13/18.
//

#ifndef MATSLISE_SE2D_H
#define MATSLISE_SE2D_H


#include <functional>
#include <vector>
#include "matslise.h"
#include "matscs.h"
#include "schrodinger.h"

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
    }

    template<int n>
    class SEnD;

    template<int n>
    struct dim {
    };

    template<>
    struct dim<1> {
        typedef std::function<double(double)> function;
        typedef ArrayXd array;
        typedef Matslise SEsolver;
    };

    template<>
    struct dim<2> {
        typedef std::function<double(double, double)> function;
        typedef ArrayXXd array;
        typedef SEnD<2> SEsolver;
    };

    template<>
    struct dim<3> {
        typedef std::function<double(double, double, double)> function;
        typedef Tensor<double, 3> array;
        typedef SEnD<3> SEsolver;
    };

    template<int n = 2>
    class SEnD {
    public:
        class Sector;

        typename dim<n>::function V;
        MatrixXd *M;
        Rectangle<n> domain;
        int sectorCount;
        typename SEnD<n>::Sector **sectors;
        ArrayXd grid[n];
        int N;
        double match;
        Options<n> options;
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

        Y<Dynamic> propagate(double E, const Y<Dynamic> &y0, double a, double b, bool use_h = true) const;

        double findEigenvalue(double Eguess, double tolerance = 1e-9, int maxIterations = 30,
                              double minTolerance = 1e-5) const;

        std::vector<typename dim<2>::array> *
        computeEigenfunction(double E, const Eigen::ArrayXd (&xs)[n]) const;

        std::vector<double> *computeEigenvaluesByIndex(int Imin, int Imax) const;
        // std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

        bool contains(double point) const {
            return point <= domain.max && point >= domain.min;
        }

        virtual ~SEnD();

    protected:
        Y<Dynamic> *computeEigenfunctionSteps(double E) const;

        MatrixXd calculateM(int k) const;

    public:
        class Sector {
        public:
            SEnD<n> *se2d;
            typename dim<n - 1>::SEsolver *matslise;
            Matscs *matscs;
            typename dim<n - 1>::array vbar;
            double min, max;

            double *eigenvalues;
            typename dim<n - 1>::array *eigenfunctions;

            Sector(SEnD<n> *se2d, double min, double max, bool backward);

            Y<Eigen::Dynamic> propagate(double E, const Y<Eigen::Dynamic> &y0, double a, double b,
                                        bool use_h = true) const;

            bool contains(double point) const {
                return point <= max && point >= min;
            }

            typename dim<n - 1>::array computeEigenfunction(int index, const ArrayXd &x) const;

            double calculateError() const;

            virtual ~Sector();

        private:
            Eigen::MatrixXd calculateDeltaV(double y) const;
        };
    };
}

#endif //MATSLISE_SE2D_H
