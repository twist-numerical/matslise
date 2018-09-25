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
#include <matslise/se2d.h>
#include <matslise/Evaluator.h>

namespace py = pybind11;
using namespace Eigen;
using namespace std;
using namespace matslise;

Y<double> packY(const tuple<double, double> &t) {
    return Y<double>({get<0>(t), get<1>(t)});
}

Y<MatrixXd> packY(const tuple<MatrixXd, MatrixXd> &t) {
    int rows = get<0>(t).rows(), cols = get<0>(t).rows();
    return Y<MatrixXd>({get<0>(t), get<1>(t)}, {MatrixXd::Zero(rows, cols), MatrixXd::Zero(rows, cols)});
}

Y<double> packY(const tuple<double, double> &t, const tuple<double, double> &dt) {
    return Y<double>({get<0>(t), get<1>(t)}, {get<0>(dt), get<1>(dt)});
}

template<typename D>
tuple<D, D> unpackY(const Y<D> &y) {
    return make_tuple(y.y[0], y.y[1]);
};

template <class R, class... A>
class PyEvaluator : public Evaluator<R, A...> {
public:
    /* Trampoline (need one for each virtual function) */
    R eval(A... x) const override {}
};

// @formatter:off
PYBIND11_MODULE(pyslise, m) {
    py::class_<Evaluator<Y<double>, double>, PyEvaluator<Y<double>, double>>(m, "EvaluatorY1")
        .def("__call__", [](Evaluator<Y<double>, double> &e, double x) -> tuple<double, double> {
            return unpackY(e(x));
        }, py::is_operator());

    py::class_<SE2D>(m, "PySE2d")
        .def(py::init<function<double(double, double)>, double,double, double,double, int, int>())
        .def("calculateError", &SE2D::calculateError)
        .def("calculateErrorMatrix", &SE2D::calculateErrorMatrix)
        .def("computeEigenfunction", &SE2D::computeEigenfunction);

    py::class_<HalfRange>(m, "PysliseHalf")
        .def(py::init<function<double(double)>, double, int>())
        .def("computeEigenvalues", [](HalfRange &m, double Emin, double Emax, const Vector2d &side) -> vector<pair<unsigned int, double>>* {
            return m.computeEigenvalues(Emin, Emax, Y<double>(side));
        })
        .def("computeEigenvaluesByIndex", [](HalfRange &m, unsigned int Imin, unsigned int Imax, const Vector2d &side) -> vector<pair<unsigned int, double>>* {
            return m.computeEigenvaluesByIndex(Imin, Imax, Y<double>(side));
        })
        .def("computeEigenfunction", [](HalfRange &m, double E, const Vector2d &side, const ArrayXd &xs)
            -> tuple<ArrayXd, ArrayXd> {
                auto ysY = m.computeEigenfunction(E, Y<double>(side), xs);
                ArrayXd ys(ysY.size());
                ArrayXd dys(ysY.size());
                for(int i = 0; i < ysY.size(); ++i) {
                    ys[i] = ysY[i].y[0];
                    dys[i] = ysY[i].y[1];
                }
                return make_tuple(ys, dys);
        });

    py::class_<Matslise>(m, "Pyslise")
        .def(py::init<function<double(double)>, double, double, int>())
        .def("propagate",
            [](Matslise &m, double E, const Vector2d &y, double a, double b) ->
                tuple<Vector2d, double> {
                    Y<double> y0;
                    double theta;
                    tie(y0, theta) = m.propagate(E, Y<double>(y), a, b);
                    return make_tuple(y0.y, theta);
                },
            py::arg("E"), py::arg("y"), py::arg("a"), py::arg("b"))
        .def("propagate",
            [](Matslise &m, double E, const Vector2d &y, const Vector2d &dy, double a, double b)
            -> tuple<Vector2d, Vector2d, double> {
                Y<double> y0;
                double theta;
                tie(y0, theta) = m.propagate(E, Y<double>(y, dy), a, b);
                return make_tuple(y0.y, y0.dy, theta);
            },
            py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
        .def("computeEigenvalues",
            [](Matslise &m, double Emin, double Emax, const Vector2d &left, const Vector2d &right)
            -> vector<pair<unsigned int, double>>* {
                return m.computeEigenvalues(Emin, Emax, Y<double>(left), Y<double>(right));
            },
            py::arg("Emin"), py::arg("Emax"), py::arg("left"), py::arg("right"))
        .def("computeEigenvaluesByIndex",
            [](Matslise &m, unsigned int Imin, unsigned int Imax, const Vector2d &left, const Vector2d &right)
            -> vector<pair<unsigned int, double>>* {
                return m.computeEigenvaluesByIndex(Imin, Imax, Y<double>(left), Y<double>(right));
            },
            py::arg("Imin"), py::arg("Imax"), py::arg("left"), py::arg("right"))
        .def("computeEigenfunction", [](Matslise &m, double E, const Vector2d &left, const Vector2d &right, const ArrayXd &xs)
             -> tuple<ArrayXd, ArrayXd> {
                auto ysY = m.computeEigenfunction(E, Y<double>(left), Y<double>(right), xs);
                ArrayXd ys(ysY.size());
                ArrayXd dys(ysY.size());
                for(int i = 0; i < ysY.size(); ++i) {
                    ys[i] = ysY[i].y[0];
                    dys[i] = ysY[i].y[1];
                }
                return make_tuple(ys, dys);
            },
            py::arg("E"), py::arg("left"), py::arg("right"), py::arg("xs"))
        .def("calculateError",
            [](Matslise &m, double E, const Vector2d &left, const Vector2d &right)
            -> tuple<double, double, double> {
                return m.calculateError(E, Y<double>(left), Y<double>(right));
            },
            py::arg("E"), py::arg("left"), py::arg("right"))
        .def("eigenfunctionCalculator",
            [](Matslise &m, double E, const Vector2d &left, const Vector2d &right) -> Evaluator<Y<double>, double>* {
                return m.eigenfunctionCalculator(E, Y<double>(left), Y<double>(right));
            },
            py::arg("E"), py::arg("left"), py::arg("right"));

    py::class_<Matscs>(m, "Pyscs")
        .def(py::init<function<MatrixXd(double)>, int, double, double, int>())
        .def("propagate",
            [](Matscs &m, double E, tuple<MatrixXd, MatrixXd> y, double a, double b) -> tuple<MatrixXd, MatrixXd> {
                return unpackY(m.propagate(E, packY(y), a, b));
        });
}