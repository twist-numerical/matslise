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
        int N;
        int gridPoints;

    public:
        SE2D(std::function<double(double, double)> V,
             double xmin, double xmax, double ymin, double ymax,
             int xSectorCount, int ySectorCount, int N = 12, int matscs_count = 5, int gridPoints=60);

        matslise::Y<MatrixXd> propagate(double E, double y, bool forward) const;

        std::pair<MatrixXd, MatrixXd> calculateErrorMatrix(double E) const;

        std::pair<double, double> calculateError(double E) const;

        std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

        std::vector<Eigen::ArrayXXd>
        computeEigenfunction(double E, const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const;

        Y<MatrixXd> *computeEigenfunctionSteps(double E) const;

        MatrixXd calculateM(int k) const;

        const se2d_util::Sector &getSector(double y) const;

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
            double *eigenfunctionsScaling;
            Eigen::ArrayXd *eigenfunctions;

            Sector(SE2D *se2d, double ymin, double ymax, int matslise_count, int matscs_count);

            matslise::Y<MatrixXd> propagate(double E, const matslise::Y<MatrixXd> &c, bool forward) const;

            matslise::Y<MatrixXd> propagate(double E, const matslise::Y<MatrixXd> &c, double y, bool forward) const;

            ArrayXd computeEigenfunction(int index, const ArrayXd &x) const;

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
