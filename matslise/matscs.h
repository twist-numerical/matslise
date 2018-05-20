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
#include <Eigen/Dense>

#define MATSCS_HMAX 10
#define MATSCS_ETA 4
#define MATSCS_LEGENDRE 6

using namespace Eigen;

namespace matscs {
    class Sector;

    class Y;
}

class Matscs {
public:
    std::function<MatrixXd(double)> V;
    int n;
    double xmin, xmax;
    int sectorCount;
    matscs::Sector **sectors;
    MatrixXd zero, one;
public:
    Matscs(std::function<MatrixXd(double)> V, int n, double xmin, double xmax, int sectorCount);

    matscs::Y propagate(double E, const matscs::Y &y, double a, double b) const;

    std::vector<matscs::Y> *computeEigenfunction(double E, std::vector<double> &x) const;

    virtual ~Matscs();
};

namespace matscs {
    class Y {
    public:
        MatrixXd y, dy;

        Y(MatrixXd y, MatrixXd dy) : y(y), dy(dy) {}

        friend std::ostream &operator<<(std::ostream &os, const Y &m) {
            return os << "(" << m.y << "," << m.dy << ")";
        }
    };

    class T {
    public:
        MatrixXd u, up, v, vp;

        T(MatrixXd u, MatrixXd up, MatrixXd v, MatrixXd vp) : u(u), up(up), v(v), vp(vp) {}

        Y operator*(const Y &y) const {
            return Y(u * y.y + v * y.dy, up * y.y + vp * y.dy);
        }

        Y operator/(const Y &y) const {
            MatrixXd d = (u * vp - up * v).inverse();
            return Y(vp * d * y.y - v * d * y.dy, -up * d * y.y + u * d * y.dy);
        }

        friend std::ostream &operator<<(std::ostream &os, const T &m) {
            return os << "((" << m.u << "," << m.v << "), (" << m.up << "," << m.vp << "))";
        }
    };

    class Sector {
    private:
        const Matscs *s;
    public:
        MatrixXd *vs;
        double xmin, xmax, h;
        Array2D<MatrixXd, MATSCS_ETA, MATSCS_HMAX> u;
        Array2D<MatrixXd, MATSCS_ETA, MATSCS_HMAX> up;
        Array2D<MatrixXd, MATSCS_ETA, MATSCS_HMAX> v;
        Array2D<MatrixXd, MATSCS_ETA, MATSCS_HMAX> vp;
        MatrixXd hu[MATSCS_ETA];
        MatrixXd hup[MATSCS_ETA];
        MatrixXd hv[MATSCS_ETA];
        MatrixXd hvp[MATSCS_ETA];
        MatrixXd D;

        Sector(const Matscs *problem, double xmin, double xmax);

        void calculateTCoeffs();

        T calculateT(double E) const;

        T calculateT(double E, double delta) const;

        virtual ~Sector();

    };
}

#endif