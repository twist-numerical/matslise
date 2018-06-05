#include <pybind11/pybind11.h>
#include <matslise.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <matscs.h>
#include <tuple>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

namespace py = pybind11;
using namespace Eigen;
using namespace std;

matslise::Y packY(tuple<double, double> &t) {
    return matslise::Y(get<0>(t), get<1>(t));
}

tuple<double, double> unpackY(matslise::Y y) {
    tuple<double, double> t(y.y, y.dy);
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
    py::class_<Matslise>(m, "Pyslise")
        .def(py::init<function<double(double)>, double, double, int>())
        .def("propagate",
            [](Matslise &m, double E, tuple<double, double> y, double a, double b) ->
                tuple<double, double, double> {
                    matslise::Y y0;
                    double theta;
                    tie(y0, theta) = m.propagate(E, packY(y), a, b);
                    return std::make_tuple(y0.y, y0.dy, theta);
                })
        .def("computeEigenfunction", [](Matslise &m, double E, vector<double> xs)
             -> tuple<vector<double>, vector<tuple<double, double>>*> {
                auto ysY = m.computeEigenfunction(E, xs);
                auto ys = new vector<tuple<double, double>>();
                for(matslise::Y y : *ysY)
                    ys->push_back(unpackY(y));
                delete ysY;
                tuple<vector<double>, vector<tuple<double, double>>*> xsys(xs, ys);
                return xsys;
            })
        .def("calculateError",
            [](Matslise &m, double E, tuple<double, double> left, tuple<double, double> right) ->
                double {
                    return m.calculateError(E, packY(left), packY(right));
                });

    py::class_<Matscs>(m, "Pyscs")
        .def(py::init<function<MatrixXd(double)>, int, double, double, int>())
        .def("propagate",
            [](Matscs &m, double E, tuple<MatrixXd, MatrixXd> y, double a, double b) -> tuple<MatrixXd, MatrixXd> {
                auto packedY = packY(y);
                return unpackY(m.propagate(E, packedY, a, b));
        });
}