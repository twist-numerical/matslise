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
#include "Array2D.h"
#include "Matrix2D.h"
#include "Evaluator.h"

#define MATSLISE_HMAX_delta 11
#define MATSLISE_ETA_delta  5
#define MATSLISE_ETA_h 7
#define MATSLISE_N 17

using namespace Eigen;

namespace matslise {

    template<typename D=double>
    class Y {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Vector2D<D> y, dy;

        Y();

        Y(Vector2D<D> y) : y(y), dy({0, 0}) {}


        Y(Vector2D<D> y, Vector2D<D> dy) : y(y), dy(dy) {}

        Y<D> operator-() const {
            return Y<D>(-y, -dy);
        }

        friend std::ostream &operator<<(std::ostream &os, const Y<D> &m) {
            return os << "(" << m.y[0] << "," << m.y[1] << ")" << "(" << m.dy[0] << "," << m.dy[1] << ")";
        }

        template<typename R>
        Y<D> operator*(const R &f) {
            return Y<D>(y * f, dy * f);
        }

        template<typename R>
        Y<D> &operator*=(const R &f) {
            y *= f;
            dy *= f;
            return *this;
        }

        friend Y<D> operator*(const D &f, const Y<D> &y) {
            return Y<D>(f * y.y, f * y.dy);
        }
    };

    template<typename D=double>
    class T {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Matrix2D<D> t;
        Matrix2D<D> dt;

        T() : t({D(1), D(0), D(0), D(1)}), dt({D(0), D(0), D(0), D(0)}) {}

        T(Matrix2D<D> t, Matrix2D<D> dt) : t(t), dt(dt) {}

        template<typename R=D>
        Y<R> operator*(Y<R> y) const {
            return Y<R>(t * y.y, t * y.dy + dt * y.y);
        }

        template<typename R=D>
        Y<R> operator/(Y<R> y) const {
            Matrix2D<D> ti, dti;
            ti = {t.d, -t.b, -t.c, t.a};
            dti = {dt.d, -dt.b, -dt.c, dt.a};
            return Y<R>(ti * y.y, ti * y.dy + dti * y.y);
        }

        friend std::ostream &operator<<(std::ostream &os, const T<D> &m) {
            return os << m.t << ", " << m.dt;
        }
    };

    namespace matslise_util {
        class Sector;

        class EigenfunctionCalculator;
    }

    class Matslise {
    public:
        std::function<double(double)> V;
        double xmin, xmax;
        int sectorCount;
        double match;
        matslise::matslise_util::Sector **sectors;

    public:
        Matslise(std::function<double(double)> V, double xmin, double xmax, int sectorCount);

        std::pair<matslise::Y<double>, double>
        propagate(double E, const matslise::Y<double> &y, double a, double b) const;

        Eigen::Array<matslise::Y<double>, Eigen::Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<double> &left, const matslise::Y<double> &right,
                             const Eigen::ArrayXd &x) const;

        std::tuple<double, double, double>
        calculateError(double E, const matslise::Y<double> &left, const matslise::Y<double> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<double> &left,
                           const matslise::Y<double> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<double> &left,
                                  const matslise::Y<double> &right) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, int Imin, int Imax,
                           const matslise::Y<double> &left,
                           const matslise::Y<double> &right) const;

        matslise::matslise_util::EigenfunctionCalculator *eigenfunctionCalculator(
                double E, const matslise::Y<double> &left, const matslise::Y<double> &right);

        virtual ~Matslise();
    };

    class HalfRange {
    public:
        const Matslise *ms;

    public:
        HalfRange(std::function<double(double)> V, double xmax, int sectorCount);

        Array<matslise::Y<double>, Dynamic, 1>
        computeEigenfunction(double E, const matslise::Y<double> &side, const ArrayXd &x) const;

        std::vector<std::pair<int, double>> *
        computeEigenvalues(double Emin, double Emax, const matslise::Y<double> &side) const;

        std::vector<std::pair<int, double>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<double> &side) const;

        virtual ~HalfRange();
    };

    namespace matslise_util {
        class Sector {
        public:
            Matslise *s;
            double *vs;
            double xmin, xmax, h;
            Array2D<Matrix2D<double>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> t_coeff;
            Matrix2D<double> t_coeff_h[MATSLISE_ETA_h];

            Sector(Matslise *problem, double xmin, double xmax);

            void calculateTCoeffs();

            T<double> calculateT(double E) const;

            T<double> calculateT(double E, double delta) const;

            Y<double> propagate(double E, const Y<double> &y0, bool forward) const;

            Y<double> propagate(double E, const Y<double> &y0, double delta) const;

            Y<double> propagate(double E, const Y<double> &y0, double delta, double &theta) const;

            double prufer(double E, double delta, const Y<double> &y0, const Y<double> &y1) const;

            virtual ~Sector();
        };

        class EigenfunctionCalculator : public Evaluator<Y<double>, double> {
        public:
            Matslise *ms;
            double E;
            std::vector<Y<double>> ys;

            EigenfunctionCalculator(Matslise *ms, double E, const Y<double> &left, const Y<double> &right);

            virtual Y<double> eval(double x) const;
        };
    };
}


#endif //SCHRODINGER_MATSLISE_H