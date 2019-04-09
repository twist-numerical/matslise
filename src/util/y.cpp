//
// Created by toon on 4/9/19.
//

#include "y.h"

using namespace matslise;
using namespace Eigen;

template<>
template<>
Y<-1, -1> T<-1>::operator/(Y<-1, -1> y) const {
    HouseholderQR<Matrix<double, -1, -1>> qr_t = t.householderQr();
    Matrix<double, -1, -1> yinv = qr_t.solve(y.y);
    return Y<-1, -1>(yinv, qr_t.solve(y.dy) - qr_t.solve(dt * yinv));
}

template<>
template<>
Y<-1, 1> T<-1>::operator/(Y<-1, 1> y) const {
    HouseholderQR<Matrix<double, -1, -1>> qr_t = t.householderQr();
    Matrix<double, -1, 1> yinv = qr_t.solve(y.y);
    return Y<-1, 1>(yinv, qr_t.solve(y.dy) - qr_t.solve(dt * yinv));
}

template<>
template<>
Y<1, 1> T<1>::operator/(Y<1, 1> y) const {
    Matrix2d tInv, dtInv;
    tInv << t(1, 1), -t(0, 1), -t(1, 0), t(0, 0);
    dtInv << dt(1, 1), -dt(0, 1), -dt(1, 0), dt(0, 0);
    return Y<1, 1>(tInv * y.y, tInv * y.dy + dtInv * y.y);
}