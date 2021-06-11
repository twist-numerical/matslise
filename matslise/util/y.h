#ifndef MATSLISE_Y_H
#define MATSLISE_Y_H

#include "eigen.h"
#include "scoped_timer.h"
#include <iostream>

namespace matslise {
    enum YDiff {
        None = 0, dX = 1, dE = 2, dXdE = 3
    };

    template<typename Scalar=double, int _dimension = 1, int _cols = _dimension>
    class Y {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        static constexpr int n4 = (_dimension == Eigen::Dynamic ? _dimension : 4 * _dimension);
        static constexpr int n2 = (_dimension == Eigen::Dynamic ? _dimension : 2 * _dimension);
        Eigen::Matrix<Scalar, n4, _cols> data;

        static Y<Scalar, _dimension, _dimension> Dirichlet(Eigen::Index N = _dimension) {
            Y<Scalar, _dimension, _dimension> y(N);
            y.block(dX) = Eigen::Matrix<Scalar, _dimension, _dimension>::Identity(N, N);
            return y;
        }

        static Y<Scalar, _dimension, _dimension> Neumann(Eigen::Index N = _dimension) {
            Y<Scalar, _dimension, _dimension> y(N);
            y.block(None) = Eigen::Matrix<Scalar, _dimension, _cols>::Identity(N, N);
            return y;
        }

        void reverse() {
            block(dX) *= -1;
            block(dXdE) *= -1;
        }

        explicit Y(Eigen::Index N = _dimension, Eigen::Index R = _dimension) {
            if (R == -1)
                R = N;
            if (N != -1) {
                data = Eigen::Matrix<Scalar, n4, _cols>::Zero(4 * N, R);
            }
        }

        explicit Y(const Eigen::Matrix<Scalar, n4, _cols> &data) : data(data) {}

        explicit Y(const Eigen::Matrix<Scalar, n2, _cols> &_y, const Eigen::Matrix<Scalar, n2, _cols> &_ydE) {
            assert(_y.rows() == _ydE.rows() && _y.cols() == _ydE.cols());
            Eigen::Index N = _y.rows() / 2;
            data = Eigen::Matrix<Scalar, n4, _cols>::Zero(4 * N, _y.cols());
            y() = _y;
            ydE() = _ydE;
        }

        Eigen::Index dimension() const {
            return data.rows() / 4;
        };

        Eigen::Index cols() const {
            return data.cols();
        };

        Y<Scalar, _dimension, _cols> operator-() const {
            return Y<Scalar, _dimension, _cols>(-data);
        }

        Eigen::Block<Eigen::Matrix<Scalar, n4, _cols>, _dimension, _cols> block(YDiff d = None) {
            return data.template block<_dimension, _cols>(((int) d) * dimension(), 0, dimension(), cols());
        }

        Eigen::Block<const Eigen::Matrix<Scalar, n4, _cols>, _dimension, _cols> block(YDiff d = None) const {
            return data.template block<_dimension, _cols>(((int) d) * dimension(), 0, dimension(), cols());
        }

        Eigen::Block<Eigen::Matrix<Scalar, n4, _cols>, n2, _cols> y() {
            return data.template topRows<n2>(2 * dimension());
        }

        Eigen::Block<const Eigen::Matrix<Scalar, n4, _cols>, n2, _cols> y() const {
            return data.template topRows<n2>(2 * dimension());
        }

        Eigen::Block<Eigen::Matrix<Scalar, n4, _cols>, n2, _cols> ydE() {
            return data.template bottomRows<n2>(2 * dimension());
        }

        Eigen::Block<const Eigen::Matrix<Scalar, n4, _cols>, n2, _cols> ydE() const {
            return data.template bottomRows<n2>(2 * dimension());
        }

        friend std::ostream &operator<<(std::ostream &os, const Y<Scalar, _dimension, _cols> &m) {
            return os << "(" << m.block(None) << "," << m.block(dX) << ")" << "(" << m.block(dE) << "," << m.block(dXdE)
                      << ")";
        }

        Y<Scalar, _dimension, _cols> operator*(const Scalar &f) const {
            return Y<Scalar, _dimension, _cols>(data * f);
        }

        Y<Scalar, _dimension, _cols> &operator*=(const Scalar &f) {
            data *= f;
            return *this;
        }

        Y<Scalar, _dimension, _cols> &operator*=(const Eigen::Matrix<Scalar, _cols, _cols> &M) {
            data *= M;
            return *this;
        }

        bool operator==(const Y &rhs) const {
            return data == rhs.data;
        }

        bool operator!=(const Y &rhs) const {
            return !(rhs == *this);
        }

        Y<Scalar, _dimension, 1> col(Eigen::Index j) const {
            return Y<Scalar, _dimension, 1>(data.col(j));
        }
    };

    template<typename Scalar, int n, int cols, int count>
    Y<Scalar, n, count> operator*(const Y<Scalar, n, cols> &y, const Eigen::Matrix<Scalar, cols, count> &M) {
        Y<Scalar, n, count> result(y.dimension(), M.cols());
        result.data = y.data * M;
        return result;
    }

    template<typename Scalar, int n, int cols, int rows>
    Y<Scalar, rows, cols> operator*(const Eigen::Matrix<Scalar, rows, n> &M, const Y<Scalar, n, cols> &y) {
        Y<Scalar, rows, cols> result;
        result.data = kroneckerProduct(Eigen::Matrix<Scalar, 4, 4>::Identity(), M) * y.data;
        return result;
    }

    template<typename Scalar, int n = 1, int n2 = (n == Eigen::Dynamic ? n : 2 * n)>
    class T {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Eigen::Matrix<Scalar, n2, n2> t;
        Eigen::Matrix<Scalar, n2, n2> dt;

        explicit T(int N = n) {
            if (N != -1) {
                t = Eigen::Matrix<Scalar, n2, n2>::Identity(2 * N, 2 * N);
                dt = Eigen::Matrix<Scalar, n2, n2>::Zero(2 * N, 2 * N);
            }
        }

        T(Eigen::Matrix<Scalar, n2, n2> t, Eigen::Matrix<Scalar, 2 * n, 2 * n> dt) : t(t), dt(dt) {}

        template<int r = 1>
        Y<Scalar, n, r> operator*(const Y<Scalar, n, r> &y) const {
            return Y<Scalar, n, r>(t * y.y(), t * y.ydE() + dt * y.y());
        }

        template<int r = 1>
        Y<Scalar, n, r> operator/(const Y<Scalar, n, r> &y) const {
            if constexpr (n == 1) {
                Eigen::Matrix<Scalar, 2, 2> tInv, dtInv;
                tInv << t(1, 1), -t(0, 1), -t(1, 0), t(0, 0);
                dtInv << dt(1, 1), -dt(0, 1), -dt(1, 0), dt(0, 0);
                return Y<Scalar, n, r>(tInv * y.y(), tInv * y.ydE() + dtInv * y.y());
            } else {
                MATSLISE_SCOPED_TIMER("T/Y dynamic");
                Eigen::HouseholderQR<Eigen::Matrix<Scalar, n2, n2>> qr_t = t.householderQr();
                Eigen::Matrix<Scalar, n2, r> yinv = qr_t.solve(y.y());
                return Y<Scalar, n, r>(yinv, qr_t.solve(y.ydE()) - qr_t.solve(dt * yinv));
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
