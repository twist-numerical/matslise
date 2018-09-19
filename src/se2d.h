//
// Created by toon on 6/13/18.
//

#ifndef MATSLISE_SE2D_H
#define MATSLISE_SE2D_H


#include <functional>
#include <tuple>
#include <vector>
#include "matslise.h"
#include "matscs.h"


#define N 12
#define GRID_POINTS 60

namespace matslise {
    namespace se2d_util {
        class Sector;

        class EigenfunctionCalculator;
    }

    class SE2D {
    public:
        std::function<double(double, double)> V;
        MatrixXd *M;
        double xmin, xmax, ymin, ymax;
        int sectorCount;
        se2d_util::Sector **sectors;
        ArrayXd xGrid;

    public:
        SE2D(std::function<double(double, double)> V,
             double xmin, double xmax, double ymin, double ymax, int xSectorCount, int ySectorCount);

        matslise::Y<MatrixXd> propagate(double E, double y, bool forward) const;

        MatrixXd calculateError(double E) const;

        std::vector<double> *computeEigenvalues(double Emin, double Emax) const;


        Eigen::ArrayXXd
        computeEigenfunction(double E, const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const;

        MatrixXd calculateM(int k) const;

        virtual ~SE2D();

    };

    namespace se2d_util {
        class Sector {
        public:
            SE2D *se2d;
            Matslise *matslise;
            Matscs *matscs;
            ArrayXd vbar;
            double ymin, ymax;

            double *eigenvalues;
            Eigen::ArrayXd *eigenfunctions;

            Sector(SE2D *se2d, double ymin, double ymax, int sectorCount);

            matslise::Y<MatrixXd> propagate(double E, const matslise::Y<MatrixXd> &c, bool forward) const;

            matslise::Y<MatrixXd> propagate(double E, const matslise::Y<MatrixXd> &c, double y, bool forward) const;

            Eigen::MatrixXd calculateDeltaV(double y) const;

            virtual ~Sector();
        };


        class EigenfunctionCalculator : public Evaluator<double, double, double> {
        public:
            SE2D *se2d;
            double E;
            std::vector<Y<VectorXd>> ys;

            EigenfunctionCalculator(SE2D *se2d, double E);

            virtual double eval(double x, double y) const;
        };
    }
}

#endif //MATSLISE_SE2D_H
