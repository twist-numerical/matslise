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
#include <Eigen/Dense>
#include "Array2D.h"

#define MATSLISE_HMAX 8
#define MATSLISE_ETA 5

using namespace Eigen;

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
    Matslise(std::function<double(double)> V, double xmin, double xmax, int sectorCount);

    std::tuple<matslise::Y, double> propagate(double E, const matslise::Y &y, double a, double b) const;

    std::vector<matslise::Y> *computeEigenfunction(double E, const matslise::Y &left, const matslise::Y &right, std::vector<double> &x) const;

    std::tuple<double, double, double>
    calculateError(double E, const matslise::Y &left, const matslise::Y &right) const;

    std::vector<std::tuple<unsigned int, double>> *computeEigenvalues(double Emin, double Emax, const matslise::Y &left, const matslise::Y &right) const;
    std::vector<std::tuple<unsigned int, double>> *computeEigenvaluesByIndex(unsigned int Imin, unsigned int Imax, const matslise::Y &left, const matslise::Y &right) const;
    std::vector<std::tuple<unsigned int, double>> *computeEigenvalues(double Emin, double Emax, unsigned int Imin, unsigned int Imax, const matslise::Y &left, const matslise::Y &right) const;

    virtual ~Matslise();
};

namespace matslise {
    class HalfRange {
    public:
        const Matslise *ms;

    public:
        HalfRange(std::function<double(double)> V, double xmax, int sectorCount);

        std::vector<matslise::Y> *computeEigenfunction(double E, const matslise::Y &side, std::vector<double> &x) const;

        std::vector<std::tuple<unsigned int, double>> *computeEigenvalues(double Emin, double Emax, const matslise::Y &side) const;
        std::vector<std::tuple<unsigned int, double>> *computeEigenvaluesByIndex(unsigned int Imin, unsigned int Imax, const matslise::Y &side) const;

        virtual ~HalfRange();
    };

    class Y {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Vector2d y, dy;

        Y() : y(Vector2d::Zero()), dy(Vector2d::Zero()) {}

        Y(Vector2d y) : y(y), dy(Vector2d::Zero()) {}

        Y(Vector2d y, Vector2d dy) : y(y), dy(dy) {}

        Y operator-() const{
            return Y(-y, -dy);
        }

        double theta() const {
            return atan(y[0] / y[1]);
        }

        friend std::ostream &operator<<(std::ostream &os, const Y &m) {
            return os << "(" << m.y[0] << "," << m.y[1] << ")" << "(" << m.dy[0] << "," << m.dy[1] << ")";
        }
    };

    class T {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Matrix2d t;
        Matrix2d dt;

        T() : t(Matrix2d::Identity()), dt(Matrix2d::Zero()) {}

        T(Matrix2d t, Matrix2d dt) : t(t), dt(dt) {}

        Y operator*(Y y) const {
            return Y(t * y.y, dt * y.y + t * y.dy);
        }

        Y operator/(Y y) const {
            Matrix2d ti = t.inverse();
            return Y(t.inverse() * y.y, -ti * ti * dt * y.y + ti * y.dy);
        }

        friend std::ostream &operator<<(std::ostream &os, const T &m) {
            return os << m.t;
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