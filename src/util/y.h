//
// Created by toon on 4/9/19.
//

#ifndef MATSLISE_Y_H
#define MATSLISE_Y_H

#include "eigen.h"
#include <iostream>

using namespace Eigen;

namespace matslise {
    template<int n = 1, int cols = n, int n2 = (n == Eigen::Dynamic ? n : 2 * n)>
    class Y {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Matrix<double, n2, cols> y, dy;

        static Y<n, n> Dirichlet(long N = n) {
            Y<n, n> y(N);
            y.y.bottomRows(y.getN()) = Matrix<double, n, cols>::Identity(N, N);
            return y;
        }

        Y(long N = n, long R = n) {
            if (R == -1)
                R = N;
            if (N != -1) {
                y = Matrix<double, n2, cols>::Zero(2 * N, R);
                dy = Matrix<double, n2, cols>::Zero(2 * N, R);
            }
        }

        Y(Matrix<double, n2, cols> y, Matrix<double, n2, cols> dy) : y(y), dy(dy) {
        }

        long getN() const {
            return y.rows() / 2;
        };

        Y<n2, cols> operator-() const {
            return Y<n2, cols>(-y, -dy);
        }

        Block<Eigen::Matrix<double, n2, cols>> getY(int derivative) {
            long N = getN();
            return y.block(derivative * N, 0, N, y.cols());
        }

        Block<const Eigen::Matrix<double, n2, cols>> getY(int derivative) const {
            long N = getN();
            return y.block(derivative * N, 0, N, y.cols());
        }

        Block<Eigen::Matrix<double, n2, cols>> getdY(int derivative) {
            long N = getN();
            return dy.block(derivative * N, 0, N, dy.cols());
        }

        Block<const Eigen::Matrix<double, n2, cols>> getdY(int derivative) const {
            long N = getN();
            return dy.block(derivative * N, 0, N, dy.cols());
        }

        friend std::ostream &operator<<(std::ostream &os, const Y<n, cols> &m) {
            return os << "(" << m.y[0] << "," << m.y[1] << ")" << "(" << m.dy[0] << "," << m.dy[1] << ")";
        }

        Y<n, cols> operator*(const double &f) const {
            return Y<n, cols>(y * f, dy * f);
        }

        Y<n, cols> &operator*=(const double &f) {
            y *= f;
            dy *= f;
            return *this;
        }

        template<int count>
        Y<n, count> operator*(const Matrix<double, cols, count> &f) const {
            long N = getN();
            Y<n, count> result(N, f.cols());
            for (int d = 0; d < 2; ++d) {
                result.getY(d) = getY(d) * f;
                result.getdY(d) = getdY(d) * f;
            }
            return result;
        }

        template<int rows>
        friend Y<rows, cols> operator*(const Matrix<double, rows, n> &f, const Y<n, cols> &y) {
            Y<rows, cols> result(f.rows());
            for (int d = 0; d < 2; ++d) {
                result.getY(d) = f * y.getY(d);
                result.getdY(d) = f * y.getdY(d);
            }
            return result;
        }

        Y<n, cols> &operator*=(const Matrix<double, cols, cols> &f) {
            for (int d = 0; d < 2; ++d) {
                getY(d) *= f;
                getdY(d) *= f;
            }
            return *this;
        }
    };

    template<int n = 1, int n2 = (n == Eigen::Dynamic ? n : 2 * n)>
    class T {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Matrix<double, n2, n2> t;
        Matrix<double, n2, n2> dt;

        T(int N = n) {
            if (N != -1) {
                t = Matrix<double, n2, n2>::Identity(2 * N, 2 * N);
                dt = Matrix<double, n2, n2>::Zero(2 * N, 2 * N);
            }
        }

        T(Matrix<double, n2, n2> t, Matrix<double, 2 * n, 2 * n> dt) : t(t), dt(dt) {}

        template<int r = 1>
        Y<n, r> operator*(Y<n, r> y) const {
            return Y<n, r>(t * y.y, t * y.dy + dt * y.y);
        }

        template<int r = 1>
        Y<n, r> operator/(Y<n, r> y) const;

        Block<Matrix<double, n2, n2>> getT(int row, int col) {
            long N = t.rows() / 2;
            return t.block(row * N, col * N, N, N);
        }

        Block<Matrix<double, n2, n2>> getdT(int row, int col) {
            long N = dt.rows() / 2;
            return dt.block(row * N, col * N, N, N);
        }

        friend std::ostream &operator<<(std::ostream &os, const T<n> &m) {
            return os << m.t << ", " << m.dt;
        }
    };
}


#endif //MATSLISE_Y_H
