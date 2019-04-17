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

        template<int n>
        class Sector;

        template<int n>
        class EigenfunctionCalculator;
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
    class SEBase {
    public:
        typename dim<n>::function V;
        MatrixXd *M;
        Rectangle<n> domain;
        int sectorCount;
        SEnD_util::Sector<n> **sectors;
        ArrayXd grid[n];
        int N;
        int match;

    public:
        SEBase(typename dim<n>::function V, const Rectangle<n> &domain, const Options<n> &options);

        std::pair<MatrixXd, MatrixXd> calculateErrorMatrix(double E) const;

        std::pair<double, double> calculateError(double E, const std::function<bool(
                std::pair<double, double>,
                std::pair<double, double>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER) const;

        std::pair<std::vector<MatrixXd>, std::vector<MatrixXd>> calculateAllSteps(double E) const;

        std::vector<std::pair<double, double>> *calculateErrors(double E) const;

        std::vector<std::pair<double, double>> *sortedErrors(double E, const std::function<bool(
                std::pair<double, double>,
                std::pair<double, double>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER) const;

        double findEigenvalue(double Eguess) const;

        // std::vector<double> *computeEigenvalues(double Emin, double Emax) const;


        virtual ~SEBase();

    protected:
        Y<Dynamic> *computeEigenfunctionSteps(double E) const;

        MatrixXd calculateM(int k) const;
    };

    template<int n>
    class SEnD : public SEBase<n> {
    public:
        using SEBase<n>::SEBase;
    };

    template<>
    class SEnD<2> : public SEBase<2> {
    public:
        using SEBase<2>::SEBase;

        std::vector<typename dim<2>::array> *
        computeEigenfunction(double E, const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const;

        std::vector<double> *computeEigenvaluesByIndex(int Imin, int Imax) const;
    };

    namespace SEnD_util {
        template<int n = 2>
        class Sector {
        public:
            SEBase<n> *se2d;
            typename dim<n - 1>::SEsolver *matslise;
            Matscs *matscs;
            typename dim<n - 1>::array vbar;
            double min, max;

            double *eigenvalues;
            typename dim<n - 1>::array *eigenfunctions;

            Sector(SEBase<n> *se2d, double min, double max, const Options<n> &options);

            matslise::Y<Eigen::Dynamic> propagate(double E, const matslise::Y<Eigen::Dynamic> &c, bool forward) const;

            matslise::Y<Eigen::Dynamic>

            propagate(double E, const matslise::Y<Eigen::Dynamic> &c, double y, bool forward) const;

            typename dim<n - 1>::array computeEigenfunction(int index, const ArrayXd &x) const;

            virtual ~Sector();

        private:
            Eigen::MatrixXd calculateDeltaV(double y) const;
        };


        template<int n = 2>
        class EigenfunctionCalculator : public Evaluator<double, double, double> {
        public:
            SEnD<n> *se2d;
            double E;
            std::vector<Y<Eigen::Dynamic, 1>> ys;

            EigenfunctionCalculator(SEnD<n> *se2d, double E);

            virtual double eval(double x, double y) const;
        };
    }
}

#endif //MATSLISE_SE2D_H
