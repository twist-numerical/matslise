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

inline Y<> make_y(const Vector2d &value) {
    return Y<>(value, {0, 0});
}

Y<> packY(const tuple<double, double> &t) {
    return make_y({get<0>(t), get<1>(t)});
}

Y<Dynamic> packY(const tuple<MatrixXd, MatrixXd> &t) {
    MatrixXd a, b;
    tie(a, b) = t;
    Y<Dynamic> result(a.rows(), a.cols());
    result.getY(0) = a;
    result.getY(1) = b;
    return result;
}

Y<> packY(const tuple<double, double> &t, const tuple<double, double> &dt) {
    return Y<>({get<0>(t), get<1>(t)}, {get<0>(dt), get<1>(dt)});
}

pair<pair<MatrixXd, MatrixXd>, pair<MatrixXd, MatrixXd>> unpackY(const Y<Dynamic> &y) {
    return make_pair(make_pair(y.getY(0), y.getY(1)), make_pair(y.getdY(0), y.getdY(1)));
};

pair<pair<double, double>, pair<double, double>> unpackY(const Y<> &y) {
    return make_pair(make_pair(y.y[0], y.y[1]), make_pair(y.dy[0], y.dy[1]));
};

template<class R, class... A>
class PyEvaluator : public Evaluator<R, A...> {
public:
    /* Trampoline (need one for each virtual function) */
    R eval(A... x) const override {}
};

// @formatter:off
PYBIND11_MODULE(pyslise, m) {
    py::class_<Evaluator<Y<>, double>, PyEvaluator<Y<>, double>>(m, "EvaluatorY1")
        .def("__call__", [](Evaluator<Y<>, double> &e, double x) -> tuple<double, double> {
            return unpackY(e(x)).first;
        }, py::is_operator());

    py::class_<SEnD_util::Sector<2>>(m, "PySE2dSector") // ?? std::unique_ptr<SEnD_util::Sector<2>,py::nodelete>
        .def_property_readonly("eigenvalues", [](SEnD_util::Sector<2> &s) -> vector<double>* {
            auto l = new vector<double>(s.se2d->N);
            for(int i = 0; i < l->size(); ++i)
                l->at(i) = s.eigenvalues[i];
            return l;
        })
        .def_property_readonly("eigenfunctions", [](SEnD_util::Sector<2> &s) -> vector<ArrayXd>* {
            auto l = new vector<ArrayXd>(s.se2d->N);
            for(int i = 0; i < l->size(); ++i)
                l->at(i) = s.eigenfunctions[i];
            return l;
        });

    py::class_<SEnD<2>>(m, "PySE2d")
        .def(py::init([](function<double(double, double)> V,
                         double xmin, double xmax, double ymin, double ymax,
                         int x_count, int y_count, int N, int in_sector_count, int grid_points){
                return std::unique_ptr<SEnD<2>>(new SEnD<2>(V, {{xmin, xmax}, ymin, ymax},
                        Options<2>()
                                .sectorCount(y_count)
                                .N(N)
                                .stepsPerSector(in_sector_count)
                                .gridPoints(grid_points)
                                .nested(Options<1>().sectorCount(x_count))
                        ));
            }), "Init SE2D",
            py::arg("V"),
            py::arg("xmin"), py::arg("xmax"), py::arg("ymin"), py::arg("ymax"),
            py::arg("x_count")=17, py::arg("y_count")=17, py::arg("N")=12, py::arg("in_sector_count")=5, py::arg("grid_points")=52)
        .def_readonly("N", &SEnD<2>::N)
        .def_property_readonly("M", [](SEnD<2> &p) -> vector<MatrixXd>* {
            auto l = new vector<MatrixXd>(p.sectorCount-1);
            for(int i = 0; i < p.sectorCount-1; ++i)
                l->at(i) = p.M[i];
            return l;
        })
        .def_property_readonly("sectors", [](SEnD<2> &p) -> vector<SEnD_util::Sector<2>*>* {
            auto l = new vector<SEnD_util::Sector<2>*>(p.sectorCount);
            for(int i = 0; i < p.sectorCount; ++i)
                l->at(i) = p.sectors[i];
            return l;
        })
        .def_property_readonly("grid", [](SEnD<2> &p) -> vector<ArrayXd>* {
            auto l = new vector<ArrayXd>(2);
            for(int i = 0; i < l->size(); ++i)
                l->at(i) = p.grid[i];
            return l;
        })
        .def("calculateError", &SEnD<2>::calculateError,
            py::arg("E"), py::arg("sorter")=SEnD_util::NEWTON_RAPHSON_SORTER)
        .def("calculateErrors", &SEnD<2>::sortedErrors,
            py::arg("E"), py::arg("sorter")=SEnD_util::NEWTON_RAPHSON_SORTER)
        .def("calculateErrorMatrix", &SEnD<2>::calculateErrorMatrix)
        .def("computeEigenfunction", &SEnD<2>::computeEigenfunction)
        .def("findEigenvalue", &SEnD<2>::findEigenvalue)
        .def("calculateAllSteps", &SEnD<2>::calculateAllSteps)
        .def_readonly_static("NEWTON_RAPHSON_SORTER", &SEnD_util::NEWTON_RAPHSON_SORTER)
        .def_readonly_static("ABS_SORTER", &SEnD_util::ABS_SORTER);

    py::class_<HalfRange>(m, "PysliseHalf")
        .def(py::init<function<double(double)>, double, int>())
        .def("computeEigenvalues", [](HalfRange &m, double Emin, double Emax, const Vector2d &side) -> vector<pair<int, double>>* {
            return m.computeEigenvalues(Emin, Emax, make_y(side));
        })
        .def("computeEigenvaluesByIndex", [](HalfRange &m, int Imin, int Imax, const Vector2d &side) -> vector<pair<int, double>>* {
            return m.computeEigenvaluesByIndex(Imin, Imax, make_y(side));
        })
        .def("computeEigenfunction", [](HalfRange &m, double E, const Vector2d &side, const ArrayXd &xs, int even)
            -> tuple<ArrayXd, ArrayXd> {
                auto ysY = m.computeEigenfunction(E, make_y(side), xs, even);
                ArrayXd ys(ysY.size());
                ArrayXd dys(ysY.size());
                for(int i = 0; i < ysY.size(); ++i) {
                    ys[i] = ysY[i].y[0];
                    dys[i] = ysY[i].y[1];
                }
                return make_tuple(ys, dys);
            },
            "Compute eigenfunction (in the values xs) for a given eigenvalue. You can indicate if the eigenfunction should be even or odd.",
            py::arg("E"), py::arg("side"), py::arg("xs"), py::arg("even")=-1);

    py::class_<Matslise>(m, "Pyslise")
        .def(py::init<function<double(double)>, double, double, int>())
        .def("propagate",
            [](Matslise &m, double E, const Vector2d &y, double a, double b) ->
                tuple<Vector2d, double> {
                    Y<> y0;
                    double theta;
                    tie(y0, theta) = m.propagate(E, make_y(y), a, b);
                    return make_tuple(y0.y, theta);
                },
            py::arg("E"), py::arg("y"), py::arg("a"), py::arg("b"))
        .def("propagate",
            [](Matslise &m, double E, const Vector2d &y, const Vector2d &dy, double a, double b)
            -> tuple<Vector2d, Vector2d, double> {
                Y<> y0;
                double theta;
                tie(y0, theta) = m.propagate(E, Y<>(y, dy), a, b);
                return make_tuple(y0.y, y0.dy, theta);
            },
            py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
        .def("computeEigenvalues",
            [](Matslise &m, double Emin, double Emax, const Vector2d &left, const Vector2d &right)
            -> vector<pair<int, double>>* {
                return m.computeEigenvalues(Emin, Emax, make_y(left), make_y(right));
            },
            py::arg("Emin"), py::arg("Emax"), py::arg("left"), py::arg("right"))
        .def("computeEigenvaluesByIndex",
            [](Matslise &m, int Imin, int Imax, const Vector2d &left, const Vector2d &right)
            -> vector<pair<int, double>>* {
                return m.computeEigenvaluesByIndex(Imin, Imax, make_y(left), make_y(right));
            },
            py::arg("Imin"), py::arg("Imax"), py::arg("left"), py::arg("right"))
        .def("computeEigenfunction", [](Matslise &m, double E, const Vector2d &left, const Vector2d &right, const ArrayXd &xs)
             -> tuple<ArrayXd, ArrayXd> {
                auto ysY = m.computeEigenfunction(E, make_y(left), make_y(right), xs);
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
                return m.calculateError(E, make_y(left), make_y(right));
            },
            py::arg("E"), py::arg("left"), py::arg("right"))
        .def("eigenfunctionCalculator",
            [](Matslise &m, double E, const Vector2d &left, const Vector2d &right) -> std::function<Y<>(double)> {
                return m.eigenfunctionCalculator(E, make_y(left), make_y(right));
            },
            py::arg("E"), py::arg("left"), py::arg("right"));

    py::class_<Matscs>(m, "Pyscs")
        .def(py::init<function<MatrixXd(double)>, int, double, double, int>())
        .def("propagate",
            [](Matscs &m, double E, tuple<MatrixXd, MatrixXd> y, double a, double b) -> pair<pair<MatrixXd, MatrixXd>, pair<MatrixXd, MatrixXd>> {
                return unpackY(m.propagate(E, packY(y), a, b));
        })
        .def("propagatePsi", &Matscs::propagatePsi);
}