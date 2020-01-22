#ifndef MATSLISE_Y_H
#define MATSLISE_Y_H

#include "eigen.h"
#include <iostream>

namespace matslise {
    template<typename Scalar=double, int n = 1, int cols = n, int n2 = (n == Eigen::Dynamic ? n : 2 * n)>
    class Y {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Eigen::Matrix<Scalar, n2, cols> y, dy;

        static Y<Scalar, n, n> Dirichlet(Eigen::Index N = n) {
            Y<Scalar, n, n> y(N);
            y.y.bottomRows(y.getN()) = Eigen::Matrix<Scalar, n, cols>::Identity(N, N);
            return y;
        }

        static Y<Scalar, n, n> Neumann(Eigen::Index N = n) {
            Y<Scalar, n, n> y(N);
            y.y.topRows(y.getN()) = Eigen::Matrix<Scalar, n, cols>::Identity(N, N);
            return y;
        }

        void reverse() {
            getY(1) *= -1;
            getdY(1) *= -1;
        }

        Y(Eigen::Index N = n, Eigen::Index R = n) {
            if (R == -1)
                R = N;
            if (N != -1) {
                y = Eigen::Matrix<Scalar, n2, cols>::Zero(2 * N, R);
                dy = Eigen::Matrix<Scalar, n2, cols>::Zero(2 * N, R);
            }
        }

        Y(const Eigen::Matrix<Scalar, n2, cols> &y, const Eigen::Matrix<Scalar, n2, cols> &dy) : y(y), dy(dy) {
        }

        Eigen::Index getN() const {
            return y.rows() / 2;
        };

        Y<Scalar, n2, cols> operator-() const {
            return Y<Scalar, n2, cols>(-y, -dy);
        }

        Eigen::Block<Eigen::Matrix<Scalar, n2, cols>> getY(int derivative) {
            Eigen::Index N = getN();
            return y.block(derivative * N, 0, N, y.cols());
        }

        Eigen::Block<const Eigen::Matrix<Scalar, n2, cols>> getY(int derivative) const {
            Eigen::Index N = getN();
            return y.block(derivative * N, 0, N, y.cols());
        }

        Eigen::Block<Eigen::Matrix<Scalar, n2, cols>> getdY(int derivative) {
            Eigen::Index N = getN();
            return dy.block(derivative * N, 0, N, dy.cols());
        }

        Eigen::Block<const Eigen::Matrix<Scalar, n2, cols>> getdY(int derivative) const {
            Eigen::Index N = getN();
            return dy.block(derivative * N, 0, N, dy.cols());
        }

        friend std::ostream &operator<<(std::ostream &os, const Y<Scalar, n, cols> &m) {
            return os << "(" << m.getY(0) << "," << m.getY(1) << ")" << "(" << m.getdY(0) << "," << m.getdY(1) << ")";
        }

        Y<Scalar, n, cols> operator*(const Scalar &f) const {
            return Y<Scalar, n, cols>(y * f, dy * f);
        }

        Y<Scalar, n, cols> &operator*=(const Scalar &f) {
            y *= f;
            dy *= f;
            return *this;
        }

        Y<Scalar, n, cols> &operator*=(const Eigen::Matrix<Scalar, cols, cols> &M) {
            y *= M;
            dy *= M;
            return *this;
        }

        bool operator==(const Y &rhs) const {
            return y == rhs.y && dy == rhs.dy;
        }

        bool operator!=(const Y &rhs) const {
            return !(rhs == *this);
        }

        Y<Scalar, n, 1> col(Eigen::Index j) const {
            return Y<Scalar, n, 1>(y.col(j), dy.col(j));
        }
    };

    template<typename Scalar, int n, int cols, int count>
    Y<Scalar, n, count> operator*(const Y<Scalar, n, cols> &y, const Eigen::Matrix<Scalar, cols, count> &M) {
        Y<Scalar, n, count> result(y.getN(), M.cols());
        result.y = y.y * M;
        result.dy = y.dy * M;
        return result;
    }

    template<typename Scalar, int n, int cols, int rows>
    Y<Scalar, rows, cols> operator*(const Eigen::Matrix<Scalar, rows, n> &M, const Y<Scalar, n, cols> &y) {
        Y<Scalar, rows, cols> result;
        result.y = kroneckerProduct(Eigen::Matrix<Scalar, 2, 2>::Identity(), M) * y.y;
        result.dy = kroneckerProduct(Eigen::Matrix<Scalar, 2, 2>::Identity(), M) * y.dy;
        return result;
    }

    template<typename Scalar, int n = 1, int n2 = (n == Eigen::Dynamic ? n : 2 * n)>
    class T {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Eigen::Matrix<Scalar, n2, n2> t;
        Eigen::Matrix<Scalar, n2, n2> dt;

        T(int N = n) {
            if (N != -1) {
                t = Eigen::Matrix<Scalar, n2, n2>::Identity(2 * N, 2 * N);
                dt = Eigen::Matrix<Scalar, n2, n2>::Zero(2 * N, 2 * N);
            }
        }

        T(Eigen::Matrix<Scalar, n2, n2> t, Eigen::Matrix<Scalar, 2 * n, 2 * n> dt) : t(t), dt(dt) {}

        template<int r = 1>
        Y<Scalar, n, r> operator*(Y<Scalar, n, r> y) const {
            return Y<Scalar, n, r>(t * y.y, t * y.dy + dt * y.y);
        }

        template<int r = 1>
        Y<Scalar, n, r> operator/(Y<Scalar, n, r> y) const {
            if (n == 1) {
                Eigen::Matrix<Scalar, 2, 2> tInv, dtInv;
                tInv << t(1, 1), -t(0, 1), -t(1, 0), t(0, 0);
                dtInv << dt(1, 1), -dt(0, 1), -dt(1, 0), dt(0, 0);
                return Y<Scalar, n, r>(tInv * y.y, tInv * y.dy + dtInv * y.y);
            } else {
                Eigen::HouseholderQR<Eigen::Matrix<Scalar, n2, n2>> qr_t = t.householderQr();
                Eigen::Matrix<Scalar, n2, r> yinv = qr_t.solve(y.y);
                return Y<Scalar, n, r>(yinv, qr_t.solve(y.dy) - qr_t.solve(dt * yinv));
            }
        }

        Eigen::Block<Eigen::Matrix<Scalar, n2, n2>> getT(int row, int col) {
            Eigen::Index N = t.rows() / 2;
            return t.block(row * N, col * N, N, N);
        }

        Eigen::Block<Eigen::Matrix<Scalar, n2, n2>> getdT(int row, int col) {
            Eigen::Index N = dt.rows() / 2;
            return dt.block(row * N, col * N, N, N);
        }

        friend std::ostream &operator<<(std::ostream &os, const T<Scalar, n> &m) {
            return os << m.t << ", " << m.dt;
        }
    };
}


#endif //MATSLISE_Y_H
