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
            .def_readonly("vbar", &SE2D<>::Sector::vbar)
            .def("calculateDeltaV", &SE2D<>::Sector::calculateDeltaV)
            .def_readonly("matslise", &SE2D<>::Sector::matslise, py::return_value_policy::reference)
            .def_readonly("matscs", &SE2D<>::Sector::matscs, py::return_value_policy::reference)
            .def_readonly("min", &SE2D<>::Sector::min)
            .def_readonly("max", &SE2D<>::Sector::max);

    py::class_<SE2D<>>(m, "PySE2d")
            .def(py::init([](function<double(double, double)> V,
                             double xmin, double xmax, double ymin, double ymax,
                             int x_count, double x_tol,
                             int y_count, double y_tol,
                             double tol,
                             int N, int in_sector_count, int grid_points) {
                     if (x_count != -1 && x_tol != -1) {
                         throw invalid_argument("Not both 'x_count' and 'x_tol' can be set.");
                     }
                     if (x_count == -1 && x_tol == -1) {
                         if (tol != -1)
                             x_tol = tol;
                         else
                             throw invalid_argument("One of 'x_count' and 'x_tol' must be set.");
                     }
                     if (y_count != -1 && y_tol != -1) {
                         throw invalid_argument("Not both 'y_count' and 'y_tol' can be set.");
                     }
                     if (y_count == -1 && y_tol == -1) {
                         if (tol != -1)
                             y_tol = tol;
                         else
                             throw invalid_argument("One of 'y_count' and 'y_tol' must be set.");
                     }
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
                 }),
                 R""""(\
In the __init__ function all needed data will be precomputed to effectively solve the given Schrödinger equation on the domain. Because of the precomputation the function V is only evaluated at the moment of initalisation. Calling other methods when the object is created will never evaluate V.

Note: steps along the y-axis are more computational expensive.

:param (float,float)->float V: the potential of the Schrödinger equation to solve
:param float xmin, xmax, ymin, ymax: the domain to work on.
:param float x_count, x_tolerance: Use only one of these. This is a guidance on how to pick the number of steps along the x-axis. With x_count pyslise will make uniform_steps. With x_tolerance the steps will be picked to try to keep the error lower than x_tolerance.
:param float y_count, y_tolerance: analogous the x_count and y_tolerance.
:param float tolerance: if none of x_count, x_tolerance, y_count or y_tolerance. x_tolerance and y_tolerance will be set to tolerance

The next set of parameters are more advanced and can be useful to tweak when the required accuracy isn't reached.

:param int N: the number of used basis functions on each sector. Defaults to 12.
:param int in_sector_count: the number of steps that will be taken per sector (in de y-direction). Defaults to 2.
:param in grid_points: the number of points that will be used to calculate the quadratures. Defaults to 60.
)"""",
                 py::arg("V"),
                 py::arg("xmin"), py::arg("xmax"), py::arg("ymin"), py::arg("ymax"),
                 py::arg("x_count") = -1, py::arg("x_tolerance") = -1,
                 py::arg("y_count") = -1, py::arg("y_tolerance") = -1,
                 py::arg("tolerance") = -1,
                 py::arg("N") = 12, py::arg("in_sector_count") = 2, py::arg("grid_points") = 60)
            .def("error", [](const SE2D<> &se2d, double const &E) -> pair<double, double> {
                return se2d.calculateError(E, &SEnD_util::ABS_SORTER<>);
            }, R""""(\
Compute the error given a guess for E. This error is the result of the requirement that the found eigenfunctions are continues. The error expresses how 'discontinues' the corresponding eigenfunction would be.

:param float E: the guessed eigenvalue.
:returns: A tuple with the computed error and the derivative of that error with respect to E.
)"""", py::arg("E"))
            .def("errors", &SE2D<>::calculateErrors, R""""(\
Just like PySE2d::calculateError(E) computes this function the discontinuity of the eigenfunction. The corresponding eigenfunction will be continuous once any of the N returned values is zero.

:param float E: the guessed eigenvalue.
:returns: A list of tuples with each of the computed errors and its derivative with respect to E.
)"""", py::arg("E"))
            .def("eigenfunction", &SE2D<>::computeEigenfunction, R""""(\
Compute all the corresponding eigenfunctions for a given eigenvalue. Most of the time this will return a singleton list. But it is possible that this eigenvalue has a higher multiplicity, so more eigenfunctions will be returned. On the other hand, when the given value for E isn't an eigenvalue then there doesn't exist an eigenfunction, so the returned list will be empty.

:param float E: the eigenvalue to compute eigenfunctions for.
:param [float] x y: the x and y values of the points to evaluate the eigenfunctions in.
:returns: a list of len(x) by len(y) grids of values. In each grid the value on position i, j is that eigenfunction evaluated in point x[i], y[j].
)"""", py::arg("E"), py::arg("x"), py::arg("y"))
            .def("eigenvalue", &SE2D<>::findEigenvalue, R""""(\
By using the algorithm of Newton-Raphson the closest eigenvalue around ``start`` will be searched. It keeps executing this algorithm until either the number of iterations is reached or the error drops below tolerance.

:param float start: the initial guess.
:param float tolerance: one of the stopping conditions. Defaults to 1e-9.
:param int iterations: the other stopping conditions. Defaults to 30.
:param float min_tolerance: if the maximum number of iterations is reached and the error is smaller than ``min_tolerance`` then the found value will be counted is eigenvalue. Defaults to 1e-5.
:returns: the eigenvalue found starting with ``guess``. Note that the found eigenvalue doesn't necessarily is the closest .
)"""", py::arg("start"), py::arg("tolerance") = 1e-9, py::arg("iterations") = 30, py::arg("min_tolerance") = 1e-5)
            .def("eigenvalues", &SE2D<>::findEigenvalues, R""""(\
This heuristics tries to find all the eigenvalues within a certain range [Emin, Emax]. Because this heuristics isn't an algortihm, it is certainly not certain that all eigenvalues are found. In short: the heuristics starts with a few initial guesses and tries to find all eigenvalues that it can 'see' from that first guess.

It is not a good idea to make the number of initial values large. This will increase computation time and, more importantly, it won't be necessarily better.

:param float Emin Emax: the start and end point of the range that will be searched.
:param int inital_values: the number of starting guesses that will be used. Defaults to 16.
:returns: a list of found eigenvalues. When one has a larger multiplicity it is repeated.
)"""", py::arg("Emin"), py::arg("Emax"), py::arg("initial_values") = 16)
            .def("__calculateErrorMatrix", &SE2D<>::calculateErrorMatrix)
            .def("__propagate",
                 [](const SE2D<> &m, double E, const MatrixXd &y, const MatrixXd &dy, double a, double b) ->
                         pair<MatrixXd, MatrixXd> {
                     Y<double, Dynamic> y0(m.N);
                     y0.getY(0) = y;
                     y0.getY(1) = dy;
                     return unpackY(m.propagate(E, y0, a, b)).first;
                 },
                 py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
            .def_readonly("__N", &SE2D<>::N, "The number of basis functions used on each sector")
            .def_property_readonly("__M", [](SE2D<> &p) -> vector<MatrixXd> * {
                auto l = new vector<MatrixXd>(static_cast<vector<MatrixXd>::size_type>(p.sectorCount - 1));
                for (int i = 0; i < p.sectorCount - 1; ++i)
                    l->at(i) = p.M[i];
                return l;
            })
            .def_property_readonly("__sectors", [](SE2D<> &p) -> vector<SE2D<>::Sector *> * {
                auto l = new vector<SE2D<>::Sector *>(static_cast<vector<SE2D<>::Sector *>::size_type>(p.sectorCount));
                for (int i = 0; i < p.sectorCount; ++i)
                    l->at(i) = p.sectors[i];
                return l;
            });

    py::class_<HalfRange<>>(m, "PySliseHalf")
            .def(py::init([](function<double(double)> V, double xmax, int steps, double tolerance) {
                if (steps != -1 && tolerance != -1)
                    throw invalid_argument("Not both 'steps' and 'tolerance' can be set.");
                if (steps == -1 && tolerance == -1)
                    throw invalid_argument("One of 'steps' and 'tolerance' must be set.");
                return new HalfRange<>(V, xmax,
                                       steps != -1 ? Matslise<>::UNIFORM(steps) : Matslise<>::AUTO(tolerance));
            }), "PySlise", py::arg("V"), py::arg("xmax"), py::arg("steps") = -1, py::arg("tolerance") = -1)
            .def("computeEigenvalues", [](HalfRange<double> &m, double Emin, double Emax,
                                          const Vector2d &side) -> vector<pair<int, double>> * {
                return m.computeEigenvalues(Emin, Emax, make_y(side));
            })
            .def("computeEigenvaluesByIndex",
                 [](HalfRange<double> &m, int Imin, int Imax, const Vector2d &side) -> vector<pair<int, double>> * {
                     return m.computeEigenvaluesByIndex(Imin, Imax, make_y(side));
                 })
            .def("computeEigenfunction",
                 [](HalfRange<double> &m, double E, const Vector2d &side, const ArrayXd &xs, int even)
                         -> tuple<ArrayXd, ArrayXd> {
                     auto ysY = m.computeEigenfunction(E, make_y(side), xs, even);
                     ArrayXd ys(ysY.size());
                     ArrayXd dys(ysY.size());
                     for (int i = 0; i < ysY.size(); ++i) {
                         ys[i] = ysY[i].y[0];
                         dys[i] = ysY[i].y[1];
                     }
                     return make_tuple(ys, dys);
                 },
                 "Compute eigenfunction (in the values xs) for a given eigenvalue. You can indicate if the eigenfunction should be even or odd.",
                 py::arg("E"), py::arg("side"), py::arg("xs"), py::arg("even") = -1);


    py::class_<Matslise<>::Sector>(m, "PySliseSector")
            .def_readonly("min", &Matslise<>::Sector::min)
            .def_readonly("max", &Matslise<>::Sector::max)
            .def_readonly("backward", &Matslise<>::Sector::backward);

    py::class_<Matslise<>>(m, "PySlise")
            .def(py::init([](function<double(double)> V, double xmin, double xmax, int steps, double tolerance) {
                     if (steps != -1 && tolerance != -1)
                         throw invalid_argument("Not both 'steps' and 'tolerance' can be set.");
                     if (steps == -1 && tolerance == -1)
                         throw invalid_argument("One of 'steps' and 'tolerance' must be set.");
                     return new Matslise<>(
                             V, xmin, xmax,
                             steps != -1 ? Matslise<>::UNIFORM(steps) : Matslise<>::AUTO(tolerance));
                 }), "PySlise", py::arg("V"), py::arg("xmin"), py::arg("xmax"), py::arg("steps") = -1,
                 py::arg("tolerance") = -1)
            .def_readonly("__sectorCount", &Matslise<>::sectorCount)
            .def_readonly("__match", &Matslise<>::match)
            .def_readonly("min", &Matslise<>::xmin)
            .def_readonly("max", &Matslise<>::xmax)
            .def("__sector", [](Matslise<> &p, int i) -> Matslise<>::Sector * {
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
            .def("computeEigenvalueError",
                 [](Matslise<> &m, double E, const Vector2d &left, const Vector2d &right)
                         -> double {
                     return m.computeEigenvalueError(E, make_y(left), make_y(right));
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"))
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
            .def("eigenfunction",
                 [](Matslise<> &m, double E, const Vector2d &left,
                    const Vector2d &right) -> std::function<pair<double, double>(double)> {
                     std::function<Y<>(double)> calculator = m.eigenfunctionCalculator(E, make_y(left), make_y(right));
                     return [calculator](double x) -> pair<double, double> {
                         Y<> y = calculator(x);
                         return make_pair(y.y[0], y.y[1]);
                     };
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"));


    py::class_<Matscs<>::Sector>(m, "PyScsSector")
            .def_readonly("min", &Matscs<>::Sector::min)
            .def_readonly("max", &Matscs<>::Sector::max)
            .def_readonly("backward", &Matscs<>::Sector::backward);
    py::class_<Matscs<>>(m, "PyScs")
            .def(py::init(
                    [](function<MatrixXd(double)> V, int N, double xmin, double xmax, int steps, double tolerance) {
                        if (steps != -1 && tolerance != -1)
                            throw invalid_argument("Not both 'steps' and 'tolerance' can be set.");
                        if (steps == -1 && tolerance == -1)
                            throw invalid_argument("One of 'steps' and 'tolerance' must be set.");
                        return new Matscs<>(
                                V, N, xmin, xmax,
                                steps != -1 ? Matscs<>::UNIFORM(steps) : Matscs<>::AUTO(tolerance));
                    }), "PyScs", py::arg("V"), py::arg("dimensions"), py::arg("xmin"), py::arg("xmax"),
                 py::arg("steps") = -1,
                 py::arg("tolerance") = -1)
            .def_readonly("__sectorCount", &Matscs<>::sectorCount)
            .def_readonly("__match", &Matscs<>::match)
            .def_readonly("min", &Matscs<>::xmin)
            .def_readonly("max", &Matscs<>::xmax)
            .def("propagate",
                 [](Matscs<> &m, double E, tuple<MatrixXd, MatrixXd> y, double a,
                    double b) -> pair<pair<MatrixXd, MatrixXd>, pair<MatrixXd, MatrixXd>> {
                     return unpackY(m.propagate(E, packY(y), a, b));
                 })
            .def("propagatePsi", &Matscs<>::propagatePsi)
            .def("__sector", [](Matscs<> &p, int i) -> Matscs<>::Sector * {
                return p.sectors[i];
            }, py::return_value_policy::reference);
}