//
// Created by toon on 6/13/18.
//

#ifndef MATSLISE_SE2D_H
#define MATSLISE_SE2D_H


#include <functional>
#include <vector>
#include "matslise.h"
#include "matscs.h"

namespace se2d {
    class Sector;
}

class SE2D {
public:
    std::function<double(double, double)> V;
    double xmin, xmax, ymin, ymax;
    int sectorCount;
    se2d::Sector **sectors;
public:
    SE2D(std::function<double(double, double)> V,
         double xmin, double xmax, double ymin, double ymax, int xSectorCount, int ySectorCount);

    std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

    virtual ~SE2D();

    ArrayXd xGrid;
};

namespace se2d {
    class Sector {
    public:
        SE2D *se2d;
        Matslise *matslise;
        Matscs *matscs;
        ArrayXd vbar;

        virtual ~Sector();

        double ymin, ymax;

        Sector(SE2D *se2d, double ymin, double ymax, int sectorCount);

        matscs::Y propagate(double delta, const matscs::Y &y);

        double *eigenvalues;
        Eigen::ArrayXd *eigenfunctions;
    };
}

#endif //MATSLISE_SE2D_H
