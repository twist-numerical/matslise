#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <tuple>

#include <Eigen/Dense>

#include <matslise/matslise.h>
#include <matslise/matscs.h>

namespace py = pybind11;
using namespace Eigen;
using namespace std;

matslise::Y packY(const tuple<double, double> &t) {
    return matslise::Y({get<0>(t), get<1>(t)});
}

matslise::Y packY(const tuple<double, double> &t, const tuple<double, double> &dt) {
    return matslise::Y({get<0>(t), get<1>(t)}, {get<0>(dt), get<1>(dt)});
}

tuple<double, double> unpackY(matslise::Y y) {
    tuple<double, double> t(y.y[0], y.y[1]);
    return t;
};

matscs::Y packY(tuple<MatrixXd, MatrixXd> &t) {
    return matscs::Y(get<0>(t), get<1>(t));
}

tuple<MatrixXd, MatrixXd> unpackY(matscs::Y y) {
    tuple<MatrixXd, MatrixXd> t(y.y, y.dy);
    return t;
};

// @formatter:off
PYBIND11_MODULE(pyslise, m) {
    py::class_<matslise::HalfRange>(m, "PysliseHalf")
        .def(py::init<function<double(double)>, double, int>())
        .def("computeEigenvalues", [](matslise::HalfRange &m, double Emin, double Emax, const Vector2d &side) -> vector<double>* {
            return m.computeEigenvalues(Emin, Emax, matslise::Y(side));
        })
        .def("computeEigenfunction", [](matslise::HalfRange &m, double E, const Vector2d &side, vector<double> xs)
            -> tuple<vector<double>, vector<Vector2d>*> {
            auto ysY = m.computeEigenfunction(E, matslise::Y(side), xs);
            auto ys = new vector<Vector2d>();
            for(matslise::Y y : *ysY)
                ys->push_back(y.y);
            delete ysY;
            return make_tuple(xs, ys);
        });

    py::class_<Matslise>(m, "Pyslise")
        .def(py::init<function<double(double)>, double, double, int>())
        .def("propagate",
            [](Matslise &m, double E, const Vector2d &y, double a, double b) ->
                tuple<Vector2d, double> {
                    matslise::Y y0;
                    double theta;
                    tie(y0, theta) = m.propagate(E, matslise::Y(y), a, b);
                    return make_tuple(y0.y, theta);
                })
        .def("propagate",
            [](Matslise &m, double E, const Vector2d &y, const Vector2d &dy, double a, double b) ->
            tuple<Vector2d, Vector2d, double> {
                matslise::Y y0;
                double theta;
                tie(y0, theta) = m.propagate(E, matslise::Y(y, dy), a, b);
                return make_tuple(y0.y, y0.dy, theta);
        })
        .def("computeEigenvalues", [](Matslise &m, double Emin, double Emax, const Vector2d &left, const Vector2d &right) -> vector<double>* {
            return m.computeEigenvalues(Emin, Emax, matslise::Y(left), matslise::Y(right));
        })
        .def("computeEigenfunction", [](Matslise &m, double E, const Vector2d &left, const Vector2d &right, vector<double> xs)
             -> tuple<vector<double>, vector<Vector2d>*> {
                auto ysY = m.computeEigenfunction(E, matslise::Y(left), matslise::Y(right), xs);
                auto ys = new vector<Vector2d>();
                for(matslise::Y y : *ysY)
                    ys->push_back(y.y);
                delete ysY;
                return make_tuple(xs, ys);
            })
        .def("calculateError",
            [](Matslise &m, double E, const Vector2d &left, const Vector2d &right) ->
                tuple<double, double, double> {
                    return m.calculateError(E, matslise::Y(left), matslise::Y(right));
                });

    py::class_<Matscs>(m, "Pyscs")
        .def(py::init<function<MatrixXd(double)>, int, double, double, int>())
        .def("propagate",
            [](Matscs &m, double E, tuple<MatrixXd, MatrixXd> y, double a, double b) -> tuple<MatrixXd, MatrixXd> {
                auto packedY = packY(y);
                return unpackY(m.propagate(E, packedY, a, b));
        });
}