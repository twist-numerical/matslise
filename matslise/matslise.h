//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_MATSLISE_H
#define SCHRODINGER_MATSLISE_H

#include <armadillo>
#include "Array2D.h"

#define HMAX 8
#define ETA 5


class Matslise {
public:
    struct Y {
        double y, dy;

        Y(double y, double dy) : y(y), dy(dy) {}

        friend std::ostream &operator<<(std::ostream &os, const Matslise::Y &m) {
            return os << "(" << m.y << "," << m.dy << ")";
        }
    };

    class T {
    public:
        double u, up, v, vp;

        T(double u, double up, double v, double vp) : u(u), up(up), v(v), vp(vp) {}

        Y operator*(Y y) {
            return Y(u * y.y + v * y.dy, up * y.y + vp * y.dy);
        }

        Y operator/(Y y) {
            double d = u * vp - up * v;
            return Y((vp * y.y - v * y.dy) / d, (-up * y.y + u * y.dy) / d);
        }

        friend std::ostream &operator<<(std::ostream &os, const Matslise::T &m) {
            return os << "((" << m.u << "," << m.v << "), (" << m.up << "," << m.vp << "))";
        }
    };

    class Sector {
    public:
        Matslise *s;
        double *vs;
        double xmin, xmax, h;
        Array2D<double, ETA, HMAX> u;
        Array2D<double, ETA, HMAX> up;
        Array2D<double, ETA, HMAX> v;
        Array2D<double, ETA, HMAX> vp;
        double hu[ETA];
        double hup[ETA];
        double hv[ETA];
        double hvp[ETA];

        Sector(Matslise *problem, double xmin, double xmax);

        void calculateTCoeffs();

        T calculateT(double E);

        T calculateT(double E, double delta);
    };

    double xmin, xmax;
    int sectorCount;
    Sector **sectors;

    double (*V)(double);

public:
    Matslise(double (*V)(double), double xmin, double xmax, int sectors);

    Y propagate(double E, Y, double a, double b);
};

#endif //SCHRODINGER_MATSLISE_H
