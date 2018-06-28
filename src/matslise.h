//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_MATSLISE_H
#define SCHRODINGER_MATSLISE_H

#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <functional>
#include <Eigen/Dense>
#include "Array2D.h"
#include "Matrix2D.h"

#define MATSLISE_HMAX 17
#define MATSLISE_ETA 9

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

    Eigen::Array<matslise::Y, Eigen::Dynamic, 1>
    computeEigenfunction(double E, const matslise::Y &left, const matslise::Y &right, const Eigen::ArrayXd &x) const;

    std::tuple<double, double, double>
    calculateError(double E, const matslise::Y &left, const matslise::Y &right) const;

    std::vector<std::tuple<unsigned int, double>> *
    computeEigenvalues(double Emin, double Emax, const matslise::Y &left, const matslise::Y &right) const;

    std::vector<std::tuple<unsigned int, double>> *
    computeEigenvaluesByIndex(unsigned int Imin, unsigned int Imax, const matslise::Y &left,
                              const matslise::Y &right) const;

    std::vector<std::tuple<unsigned int, double>> *
    computeEigenvalues(double Emin, double Emax, unsigned int Imin, unsigned int Imax, const matslise::Y &left,
                       const matslise::Y &right) const;

    virtual ~Matslise();
};

namespace matslise {
    class HalfRange {
    public:
        const Matslise *ms;

    public:
        HalfRange(std::function<double(double)> V, double xmax, int sectorCount);

        Array<matslise::Y, Dynamic, 1> computeEigenfunction(double E, const matslise::Y &side, const ArrayXd &x) const;

        std::vector<std::tuple<unsigned int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y &side) const;

        std::vector<std::tuple<unsigned int, double>> *
        computeEigenvaluesByIndex(unsigned int Imin, unsigned int Imax, const matslise::Y &side) const;

        virtual ~HalfRange();
    };

    class Y {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Vector2D<> y, dy;

        Y() : y({0,0}), dy({0,0}) {}

        Y(Vector2D<> y) : y(y), dy({0,0}) {}

        Y(Vector2D<> y, Vector2D<> dy) : y(y), dy(dy) {}

        Y operator-() const {
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
        Matrix2D<> t;
        Matrix2D<> dt;

        T() : t({1,0,0,1}), dt({0,0,0,0}) {}

        T(Matrix2D<> t, Matrix2D<> dt) : t(t), dt(dt) {}

        Y operator*(Y y) const {
            return Y(t * y.y, t * y.dy + dt * y.y);
        }

        Y operator/(Y y) const {
            Matrix2D<> ti, dti;
            ti = {t.d, -t.b, -t.c, t.a};
            dti = {dt.d, -dt.b, -dt.c, dt.a};
            return Y(ti * y.y, ti * y.dy + dti * y.y);
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
        Array2D<Matrix2D<>, MATSLISE_ETA, MATSLISE_HMAX> t_coeff;
        Matrix2D<> t_coeff_h[MATSLISE_ETA];

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