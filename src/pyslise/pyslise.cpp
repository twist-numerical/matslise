#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <tuple>

#include "../matslise.h"

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

Y<double, Dynamic> packY(const tuple<MatrixXd, MatrixXd> &t) {
    MatrixXd a, b;
    tie(a, b) = t;
    Y<double, Dynamic> result(a.rows(), a.cols());
    result.getY(0) = a;
    result.getY(1) = b;
    return result;
}

Y<> packY(const tuple<double, double> &t, const tuple<double, double> &dt) {
    return Y<>({get<0>(t), get<1>(t)}, {get<0>(dt), get<1>(dt)});
}

pair<pair<MatrixXd, MatrixXd>, pair<MatrixXd, MatrixXd>> unpackY(const Y<double, Dynamic> &y) {
    return make_pair(make_pair(y.getY(0), y.getY(1)), make_pair(y.getdY(0), y.getdY(1)));
}

pair<pair<double, double>, pair<double, double>> unpackY(const Y<> &y) {
    return make_pair(make_pair(y.y[0], y.y[1]), make_pair(y.dy[0], y.dy[1]));
}

// @formatter:off
PYBIND11_MODULE(pyslise, m) {
    py::class_<SE2D<>::Sector>(m, "PySE2dSector") // ?? std::unique_ptr<SEnD_util::Sector<2>,py::nodelete>
            .def_property_readonly("eigenvalues", [](SE2D<>::Sector &s) -> vector<double> * {
                auto l = new vector<double>(s.se2d->N);
                for (unsigned int i = 0; i < l->size(); ++i)
                    l->at(i) = s.eigenvalues[i];
                return l;
            })
            .def_property_readonly("eigenfunctions", [](SE2D<>::Sector &s) -> vector<ArrayXd> * {
                auto l = new vector<ArrayXd>(s.se2d->N);
                for (unsigned int i = 0; i < l->size(); ++i)
                    l->at(i) = s.eigenfunctions[i];
                return l;
            })
            .def_readonly("matslise", &SE2D<>::Sector::matslise, py::return_value_policy::reference)
            .def_readonly("matscs", &SE2D<>::Sector::matscs, py::return_value_policy::reference)
            .def_readonly("min", &SE2D<>::Sector::min)
            .def_readonly("max", &SE2D<>::Sector::max);

    py::class_<SE2D<>>(m, "PySE2d")
            .def(py::init([](function<double(double, double)> V,
                             double xmin, double xmax, double ymin, double ymax,
                             int x_count, double x_tol,
                             int y_count, double y_tol,
                             int N, int in_sector_count, int grid_points) {
                     if (x_count != -1 && x_tol != -1)
                         throw invalid_argument("Not both 'x_count' and 'x_tol' can be set.");
                     if (x_count == -1 && x_tol == -1)
                         throw invalid_argument("One of 'x_count' and 'x_tol' must be set.");
                     if (y_count != -1 && y_tol != -1)
                         throw invalid_argument("Not both 'y_count' and 'y_tol' can be set.");
                     if (y_count == -1 && y_tol == -1)
                         throw invalid_argument("One of 'y_count' and 'y_tol' must be set.");
                     Options1<> o1;
                     if (x_count != -1)
                         o1.sectorCount(x_count);
                     else
                         o1.tolerance(x_tol);
                     Options2<> o2;
                     o2.N(N)
                             .stepsPerSector(in_sector_count)
                             .gridPoints(grid_points)
                             .nested(o1);
                     if (y_count != -1)
                         o2.sectorCount(y_count);
                     else
                         o2.tolerance(y_tol);
                     return std::unique_ptr<SE2D<>>(new SE2D<>(V, {{xmin, xmax}, ymin, ymax}, o2));
                 }), "Init SE2D",
                 py::arg("V"),
                 py::arg("xmin"), py::arg("xmax"), py::arg("ymin"), py::arg("ymax"),
                 py::arg("x_count") = -1, py::arg("x_tol") = -1,
                 py::arg("y_count") = -1, py::arg("y_tol") = -1,
                 py::arg("N") = 12, py::arg("in_sector_count") = 5, py::arg("grid_points") = 52)
            .def_readonly("N", &SE2D<>::N)
            .def_property_readonly("M", [](SE2D<> &p) -> vector<MatrixXd> * {
                auto l = new vector<MatrixXd>(p.sectorCount - 1);
                for (int i = 0; i < p.sectorCount - 1; ++i)
                    l->at(i) = p.M[i];
                return l;
            })
            .def_property_readonly("sectors", [](SE2D<> &p) -> vector<SE2D<>::Sector *> * {
                auto l = new vector<SE2D<>::Sector *>(p.sectorCount);
                for (int i = 0; i < p.sectorCount; ++i)
                    l->at(i) = p.sectors[i];
                return l;
            })
            .def_property_readonly("grid", [](SE2D<> &p) -> vector<ArrayXd> * {
                auto l = new vector<ArrayXd>(2);
                for (unsigned int i = 0; i < l->size(); ++i)
                    l->at(i) = p.grid[i];
                return l;
            })
            .def("calculateError", &SE2D<>::calculateError,
                 py::arg("E"),
                 py::arg("sorter") = static_cast<function<bool(pair<double, double>, pair<double, double>)>>(
                         &SEnD_util::NEWTON_RAPHSON_SORTER<double>))
            .def("calculateErrors", &SE2D<>::sortedErrors,
                 py::arg("E"),
                 py::arg("sorter") = static_cast<function<bool(pair<double, double>, pair<double, double>)>>(
                         &SEnD_util::NEWTON_RAPHSON_SORTER<double>))
            .def("calculateErrorMatrix", &SE2D<>::calculateErrorMatrix)
            .def("computeEigenfunction",
                 [](SE2D<> &m, double E, const ArrayXd &x, const ArrayXd &y) {
                     return m.computeEigenfunction(E, x, y);
                 })
            .def("findEigenvalue", &SE2D<>::findEigenvalue, py::arg("start"), py::arg("tolerance") = 1e-9,
                 py::arg("maxIterations") = 30, py::arg("minTolerance") = 1e-5)
            .def("propagate",
                 [](const SE2D<> &m, double E, const MatrixXd &y, const MatrixXd &dy, double a, double b) ->
                         pair<MatrixXd, MatrixXd> {
                     Y<double, Dynamic> y0(m.N);
                     y0.getY(0) = y;
                     y0.getY(1) = dy;
                     return unpackY(m.propagate(E, y0, a, b)).first;
                 },
                 py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
            .def_static("NEWTON_RAPHSON_SORTER", &SEnD_util::NEWTON_RAPHSON_SORTER<>)
            .def_static("ABS_SORTER", &SEnD_util::ABS_SORTER<>);
/*
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
*/

    py::class_<Matslise<>::Sector>(m, "PysliseSector")
            .def_readonly("min", &Matslise<>::Sector::min)
            .def_readonly("max", &Matslise<>::Sector::max)
            .def_readonly("backward", &Matslise<>::Sector::backward);

    py::class_<Matslise<>>(m, "Pyslise")
            .def(py::init([](function<double(double)> V, double xmin, double xmax, int steps, double tolerance) {
                     if (steps != -1 && tolerance != -1)
                         throw invalid_argument("Not both 'steps' and 'tolerance' can be set.");
                     if (steps == -1 && tolerance == -1)
                         throw invalid_argument("One of 'steps' and 'tolerance' must be set.");
                     return new Matslise<>(
                             V, xmin, xmax,
                             steps != -1 ? Matslise<>::UNIFORM(steps) : Matslise<>::AUTO(tolerance));
                 }), "Pyslise", py::arg("V"), py::arg("xmin"), py::arg("xmax"), py::arg("steps") = -1,
                 py::arg("tolerance") = -1)
            .def_readonly("sectorCount", &Matslise<>::sectorCount)
            .def_readonly("match", &Matslise<>::match)
            .def_readonly("min", &Matslise<>::xmin)
            .def_readonly("max", &Matslise<>::xmax)
            .def("sector", [](Matslise<> &p, int i) -> Matslise<>::Sector * {
                return p.sectors[i];
            }, py::return_value_policy::reference)
            .def("propagate",
                 [](Matslise<> &m, double E, const Vector2d &y, double a, double b) ->
                         tuple<Vector2d, double> {
                     Y<> y0;
                     double theta;
                     tie(y0, theta) = m.propagate(E, make_y(y), a, b);
                     return make_tuple(y0.y, theta);
                 },
                 py::arg("E"), py::arg("y"), py::arg("a"), py::arg("b"))
            .def("propagate",
                 [](Matslise<> &m, double E, const Vector2d &y, const Vector2d &dy, double a, double b)
                         -> tuple<Vector2d, Vector2d, double> {
                     Y<> y0;
                     double theta;
                     tie(y0, theta) = m.propagate(E, Y<>(y, dy), a, b);
                     return make_tuple(y0.y, y0.dy, theta);
                 },
                 py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
            .def("computeEigenvalues",
                 [](Matslise<> &m, double Emin, double Emax, const Vector2d &left, const Vector2d &right)
                         -> vector<pair<int, double>> * {
                     return m.computeEigenvalues(Emin, Emax, make_y(left), make_y(right));
                 },
                 py::arg("Emin"), py::arg("Emax"), py::arg("left"), py::arg("right"))
            .def("computeEigenvaluesByIndex",
                 [](Matslise<> &m, int Imin, int Imax, const Vector2d &left, const Vector2d &right)
                         -> vector<pair<int, double>> * {
                     return m.computeEigenvaluesByIndex(Imin, Imax, make_y(left), make_y(right));
                 },
                 py::arg("Imin"), py::arg("Imax"), py::arg("left"), py::arg("right"))
            .def("computeEigenfunction",
                 [](Matslise<> &m, double E, const Vector2d &left, const Vector2d &right, const ArrayXd &xs)
                         -> tuple<ArrayXd, ArrayXd> {
                     auto ysY = m.computeEigenfunction(E, make_y(left), make_y(right), xs);
                     ArrayXd ys(ysY.size());
                     ArrayXd dys(ysY.size());
                     for (int i = 0; i < ysY.size(); ++i) {
                         ys[i] = ysY[i].y[0];
                         dys[i] = ysY[i].y[1];
                     }
                     return make_tuple(ys, dys);
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"), py::arg("xs"))
            .def("calculateError",
                 [](Matslise<> &m, double E, const Vector2d &left, const Vector2d &right)
                         -> tuple<double, double, double> {
                     return m.calculateError(E, make_y(left), make_y(right));
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"))
            .def("eigenfunctionCalculator",
                 [](Matslise<> &m, double E, const Vector2d &left,
                    const Vector2d &right) -> std::function<pair<double, double>(double)> {
                     std::function<Y<>(double)> calculator = m.eigenfunctionCalculator(E, make_y(left), make_y(right));
                     return [calculator](double x) -> pair<double, double> {
                         Y<> y = calculator(x);
                         return make_pair(y.y[0], y.y[1]);
                     };
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"));


    py::class_<Matscs<>::Sector>(m, "PyscsSector")
            .def_readonly("min", &Matscs<>::Sector::min)
            .def_readonly("max", &Matscs<>::Sector::max)
            .def_readonly("backward", &Matscs<>::Sector::backward);
    py::class_<Matscs<>>(m, "Pyscs")
            .def(py::init(
                    [](function<MatrixXd(double)> V, int N, double xmin, double xmax, int steps, double tolerance) {
                        if (steps != -1 && tolerance != -1)
                            throw invalid_argument("Not both 'steps' and 'tolerance' can be set.");
                        if (steps == -1 && tolerance == -1)
                            throw invalid_argument("One of 'steps' and 'tolerance' must be set.");
                        return new Matscs<>(
                                V, N, xmin, xmax,
                                steps != -1 ? Matscs<>::UNIFORM(steps) : Matscs<>::AUTO(tolerance));
                    }), "Pyslise", py::arg("V"), py::arg("dimensions"), py::arg("xmin"), py::arg("xmax"),
                 py::arg("steps") = -1,
                 py::arg("tolerance") = -1)
            .def_readonly("sectorCount", &Matscs<>::sectorCount)
            .def_readonly("match", &Matscs<>::match)
            .def_readonly("min", &Matscs<>::xmin)
            .def_readonly("max", &Matscs<>::xmax)
            .def("propagate",
                 [](Matscs<> &m, double E, tuple<MatrixXd, MatrixXd> y, double a,
                    double b) -> pair<pair<MatrixXd, MatrixXd>, pair<MatrixXd, MatrixXd>> {
                     return unpackY(m.propagate(E, packY(y), a, b));
                 })
            .def("propagatePsi", &Matscs<>::propagatePsi)
            .def("sector", [](Matscs<> &p, int i) -> Matscs<>::Sector * {
                return p.sectors[i];
            }, py::return_value_policy::reference);
}