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

#define MATSCS_HMAX 5
#define MATSCS_ETA 3

using namespace Eigen;

namespace matslise {
    namespace matscs_sector {
        class Sector;
    }

    class Matscs {
    public:
        std::function<MatrixXd(double)> V;
        int n;
        double xmin, xmax;
        int sectorCount;
        matslise::matscs_sector::Sector **sectors;
        MatrixXd zero, one;
    public:
        Matscs(std::function<MatrixXd(double)> V, int n, double xmin, double xmax, int sectorCount);

        matslise::Y<MatrixXd> propagate(double E, const matslise::Y<MatrixXd> &y, double a, double b) const;

        std::vector<matslise::Y<MatrixXd>> *computeEigenfunction(double E, std::vector<double> &x) const;

        virtual ~Matscs();
    };

    namespace matscs_sector {
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

            virtual ~Sector();

        };
    }
}

#endif