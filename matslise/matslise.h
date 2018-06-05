//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_MATSLISE_H
#define SCHRODINGER_MATSLISE_H

#include <ostream>
#include <array>
#include <vector>
#include <tuple>
#include <functional>
#include "Array2D.h"

#define MATSLISE_HMAX 8
#define MATSLISE_ETA 5

namespace matslise {
    class Sector;

    class Y;
}

class Matslise {
public:
    std::function<double(double)> V;
    double xmin, xmax;
    int sectorCount;
    double match;
    matslise::Sector **sectors;

public:
    Matslise(std::function<double(double)>, double xmin, double xmax, int sectorCount);

    std::tuple<matslise::Y, double> propagate(double E, const matslise::Y &y, double a, double b) const;

    std::vector<matslise::Y> *computeEigenfunction(double E, std::vector<double> &x) const;

    double calculateError(double E, const matslise::Y &left, const matslise::Y &right) const;

    virtual ~Matslise();
};

namespace matslise {
    class Y {
    public:
        double y, dy;

        Y() : y(0), dy(0) {}

        Y(double y, double dy) : y(y), dy(dy) {}

        double theta() const {
            return atan(y/dy);
        }

        friend std::ostream &operator<<(std::ostream &os, const Y &m) {
            return os << "(" << m.y << "," << m.dy << ")";
        }
    };

    class T {
    public:
        double u, up, v, vp;

        T(double u, double up, double v, double vp) : u(u), up(up), v(v), vp(vp) {}

        Y operator*(Y y) const {
            return Y(u * y.y + v * y.dy, up * y.y + vp * y.dy);
        }

        Y operator/(Y y) const {
            double d = u * vp - up * v;
            return Y((vp * y.y - v * y.dy) / d, (-up * y.y + u * y.dy) / d);
        }

        friend std::ostream &operator<<(std::ostream &os, const T &m) {
            return os << "((" << m.u << "," << m.v << "), (" << m.up << "," << m.vp << "))";
        }
    };

    class Sector {
    public:
        Matslise *s;
        double *vs;
        double xmin, xmax, h;
        Array2D<double, MATSLISE_ETA, MATSLISE_HMAX> u;
        Array2D<double, MATSLISE_ETA, MATSLISE_HMAX> up;
        Array2D<double, MATSLISE_ETA, MATSLISE_HMAX> v;
        Array2D<double, MATSLISE_ETA, MATSLISE_HMAX> vp;
        double hu[MATSLISE_ETA];
        double hup[MATSLISE_ETA];
        double hv[MATSLISE_ETA];
        double hvp[MATSLISE_ETA];

        Sector(Matslise *problem, double xmin, double xmax);

        void calculateTCoeffs();

        T calculateT(double E) const;

        T calculateT(double E, double delta) const;

        Y propagate(double E, const Y &y0, double delta, double &theta) const;
        double prufer(double E, double delta, const Y &y0, const Y &y1) const;
        virtual ~Sector();
    };
}


#endif //SCHRODINGER_MATSLISE_H