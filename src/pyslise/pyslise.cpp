#include <utility>


#include "module.h"

void pyslise(py::module &m) {
    py::class_<AbstractMatslise<double>>(m, "AbstractPyslise")
            .def("eigenvalues",
                 [](AbstractMatslise<double> &m, double Emin, double Emax, const Vector2d &left,
                    const optional<Vector2d> &_right)
                         -> vector<pair<int, double>> {
                     const Vector2d &right = _right ? *_right : left;
                     return m.eigenvalues(Emin, Emax, make_y(left), make_y(right));
                 }, R""""(\
Calculate the eigenvalues in an interval [Emin; Emax]. The boundary conditions have to be specified.

:param float Emin, Emax: the ends of the search interval.
:param (float,float) left, right: the boundary conditions. The corresponding eigenfunctions will have as left value a scaling of left[0] and derivative left[1].

:returns: a list of tuples. Each tuples contains the index and the eigenvalue with that index.
)"""",
                 py::arg("Emin"), py::arg("Emax"), py::arg("left"), py::arg("right") = optional<Vector2d>())
            .def("eigenvaluesByIndex",
                 [](AbstractMatslise<double> &m, int Imin, int Imax, const Vector2d &left,
                    const optional<Vector2d> &_right) -> vector<pair<int, double>> {
                     const Vector2d &right = _right ? *_right : left;
                     return m.eigenvaluesByIndex(Imin, Imax, make_y(left), make_y(right));
                 }, R""""(\
Calculate all eigenvalues with index between Imin and Imax. The first eigenvalue has index 0. Imin inclusive, Imax exclusive.

:param int Imin: the first eigenvalue to find, by index.
:param int Imax: only the first Imax eigenvalues will be considered.
:param (float,float) left, right: the boundary conditions.

:returns: a list of tuples. Each tuples contains the index and the eigenvalue with that index.
)"""",
                 py::arg("Imin"), py::arg("Imax"), py::arg("left"), py::arg("right") = optional<Vector2d>())
            .def("eigenvalueError",
                 [](AbstractMatslise<double> &m, double E, const Vector2d &left, const optional<Vector2d> &_right,
                    int index)
                         -> double {
                     const Vector2d &right = _right ? *_right : left;
                     return m.eigenvalueError(E, make_y(left), make_y(right), index);
                 }, R""""(\
Calculate the error for a given eigenvalue. It will use a less accurate method to estimate another (worse) guess for that eigenvalue. The true error on the given eigenvalue will be less than the value returned by this method.

:param float E: the eigenvalue to calculate the error for.
:param (float,float) left, right: the boundary conditions.

:returns: the error.
)"""",
                 py::arg("E"), py::arg("left"), py::arg("right") = optional<Vector2d>(), py::arg("index") = -1)
            .def("eigenfunction",
                 [](AbstractMatslise<double> &m, double E, const ArrayXd &xs, const Vector2d &left,
                    const optional<Vector2d> &_right, int index)
                         -> tuple<ArrayXd, ArrayXd> {
                     const Vector2d &right = _right ? *_right : left;
                     auto ysY = m.eigenfunction(E, make_y(left), make_y(right), index)(xs);
                     ArrayXd ys(ysY.size());
                     ArrayXd dys(ysY.size());
                     for (Eigen::Index i = 0; i < ysY.size(); ++i) {
                         ys[i] = ysY[i].y[0];
                         dys[i] = ysY[i].y[1];
                     }
                     return make_tuple(ys, dys);
                 }, R""""(\
Calculate the eigenfunction corresponding to the eigenvalue E in the points xs.

:param float E: the eigenvalue.
:param (float,float) left, right: the boundary conditions.
:param xs: the points to calculate the eigenfunction for.

:returns: a pair of lists which each a length of len(xs). The first list contains the values of the eigenfunction in the points xs. The second contains the derivative of the eigenfunction in those points.
)"""",
                 py::arg("E"), py::arg("xs"), py::arg("left"), py::arg("right") = optional<Vector2d>(),
                 py::arg("index") = -1)
            .def("eigenfunctionCalculator",
                 [](AbstractMatslise<double> &m, double E, const Vector2d &left,
                    const optional<Vector2d> &_right, int index) -> function<pair<double, double>(double)> {
                     const Vector2d &right = _right ? *_right : left;
                     function<Y<>(double)> calculator = m.eigenfunction(E, make_y(left), make_y(right), index);
                     return [calculator](double x) -> pair<double, double> {
                         Y<> y = calculator(x);
                         return make_pair(y.y[0], y.y[1]);
                     };
                 }, R""""(\
Returns the eigenfunction corresponding to the eigenvalue E as a python function. The returned function can be evaluated in all the points in the domain.

:param float E: the eigenvalue.
:param (float,float) left, [right]: the boundary conditions.
:param int [index]: the index of the eigenvalue

:returns: a function that takes a value and returns a tuple with the eigenfunction and derivative in that value.
)"""",
                 py::arg("E"), py::arg("left"), py::arg("right") = optional<Vector2d>(), py::arg("index") = -1)
            .def_property_readonly("domain", [](const AbstractMatslise<double> &matslise) {
                return std::pair<double, double>(matslise.domain.min, matslise.domain.max);
            });


    py::class_<Matslise<>, AbstractMatslise<double>>(m, "Pyslise", R""""(\
>>> from math import pi, cos
>>> p = Pyslise(lambda x: 2*cos(2*x), 0, pi, 1e-8)
>>> i, E = p.eigenvaluesByIndex(0, 1, (0,1))[0]
>>> abs(E - -0.1102488054219) < 1e-6
True
)"""")
            .def(py::init([](const function<double(double)> &V, double min, double max, double tolerance) {
                return new Matslise<>(V, min, max, tolerance);
            }), R""""(\
In the __init__ function all needed data will be precomputed to effectively solve the Schrödinger equation with given potential on the interval [min; max]. Because of the precomputation the function V is only evaluated at the moment of initalisation. Calling other methods after the object is created never will evaluate V.

Note: only one of steps and tolerance have to be set.

:param (float)->float V: the potential.
:param float min, max: the ends of the domain.
:param int steps: the number of steps to take.
:param int tolerance: automatically choose steps with at least the given accuracy.
)"""", py::arg("V"), py::arg("min"), py::arg("max"), py::arg("tolerance") = 1e-8)
            .def("propagate",
                 [](Matslise<> &m, double E, const Vector2d &y, double a, double b) ->
                         tuple<Vector2d, double> {
                     Y<> y0;
                     double theta;
                     tie(y0, theta) = m.propagate(E, make_y(y), a, b);
                     return make_tuple(y0.y, theta);
                 }, R""""(\
For a given E and initial condition in point a, propagate the solution of the Schrödinger equation to the point b.

:param float E: the fixed eigenvalue to use. This value doesn't have to be a true eigenvalue.
:param (float,float) y: the initial condition in point a.
:param float a: the start of the propagation.
:param float b: the point in which the solution is sought.

:returns: a tuple of length two. The first element is a tuple with the value and the derivative in the point b. The second element is the angle found by the Prüfer transformation.
)"""",
                 py::arg("E"), py::arg("y"), py::arg("a"), py::arg("b"))
            .def("propagate",
                 [](Matslise<> &m, double E, const Vector2d &y, const Vector2d &dy, double a, double b)
                         -> tuple<Vector2d, Vector2d, double> {
                     Y<> y0;
                     double theta;
                     tie(y0, theta) = m.propagate(E, Y<>(y, dy), a, b);
                     return make_tuple(y0.y, y0.dy, theta);
                 }, R""""(\
For a given E and initial condition in point a, propagate the solution of the Schrödinger equation to the point b. In this method the derivative with respect to E is also calculated.

:param float E: the fixed eigenvalue to use. This value doesn't have to be a true eigenvalue.
:param (float,float) y: the initial condition in point a.
:param (float,float) dy: the initial condition derived to E in point a.
:param float a: the start of the propagation.
:param float b: the point in which the solution is sought.

:returns: a tuple of length three. The first element is a tuple with the value and the derivative in the point b. The second element contains the derivative to E of the first element. The third element is the angle found by the Prüfer transformation.)"""",
                 py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
            .def("__error",
                 [](Matslise<> &m, double E, const Vector2d &left, const Vector2d &right)
                         -> tuple<double, double, double> {
                     return m.matchingError(E, make_y(left), make_y(right));
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"))
            .def_readonly("__sectorCount", &Matslise<>::sectorCount)
            .def_readonly("__matchIndex", &Matslise<>::matchIndex)
            .def("__sector", [](Matslise<> &p, int i) -> Matslise<>::Sector * {
                return p.sectors[i];
            }, py::return_value_policy::reference);

    py::class_<MatsliseHalf<>, AbstractMatslise<double>>(m, "PysliseHalf")
            .def(py::init([](const function<double(double)> &V, double xmax, double tolerance) {
                return new MatsliseHalf<>(V, xmax, tolerance);
            }), R""""(\
In the __init__ function all needed data will be precomputed to effectively solve the Schrödinger equation with given potential on the interval [min; max]. Because of the precomputation the function V is only evaluated at the moment of initalisation. Calling other methods after the object is created never will evaluate V.

Note: only one of steps and tolerance have to be set.

:param (float)->float V: the potential.
:param float min, max: the ends of the domain.
:param int steps: the number of steps to take.
:param int tolerance: automatically choose steps with at least the given accuracy.
)"""", py::arg("V"), py::arg("xmax"), py::arg("tolerance") = 1e-8);

    py::class_<Matslise<>::Sector>(m, "PysliseSector")
            .def_readonly("min", &Matslise<>::Sector::min)
            .def_readonly("max", &Matslise<>::Sector::max)
            .def_property_readonly("forward", [](const Matslise<>::Sector &s) -> bool {
                return s.direction == Direction::forward;
            })
            .def_property_readonly("v", [](const Matslise<>::Sector &s) -> vector<double> {
                vector<double> r(MATSLISE_N);
                for (int i = 0; i < MATSLISE_N; ++i)
                    r[i] = s.vs[i];
                return r;
            });

}