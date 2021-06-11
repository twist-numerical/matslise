#ifndef MATSLISE_PYSLISE
#define MATSLISE_PYSLISE

#include "../matslise.h"
#include "../matscs.h"
#include "../matslise2d.h"
#include "../matslise3d.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include <tuple>

namespace py = pybind11;
using namespace Eigen;
using namespace std;
using namespace matslise;


inline Y<> make_y(const Vector2d &value) {
    return Y<>(value, {0, 0});
}

inline Y<> packY(const tuple<double, double> &t) {
    return make_y({get<0>(t), get<1>(t)});
}

inline Y<double, Dynamic> packY(const tuple<MatrixXd, MatrixXd> &t) {
    MatrixXd a, b;
    tie(a, b) = t;
    Y<double, Dynamic> result(a.rows(), a.cols());
    result.block() = a;
    result.block(dX) = b;
    return result;
}

inline Y<> packY(const tuple<double, double> &t, const tuple<double, double> &dt) {
    return Y<>({get<0>(t), get<1>(t)},
               {get<0>(dt), get<1>(dt)});
}

inline pair<pair<MatrixXd, MatrixXd>, pair<MatrixXd, MatrixXd>> unpackY(const Y<double, Dynamic> &y) {
    return make_pair(make_pair(y.block(), y.block(dX)), make_pair(y.block(dE), y.block(dXdE)));
}

inline pair<pair<double, double>, pair<double, double>> unpackY(const Y<> &y) {
    return make_pair(make_pair(y.y()[0], y.y()[1]), make_pair(y.ydE()[0], y.ydE()[1]));
}

#endif