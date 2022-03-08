#include <utility>

#include "../liouville.h"
#include "module.h"

template<typename Scalar>
class PeriodicEigenfunctionWrapper : public PeriodicMatslise<Scalar>::Eigenfunction {
public:
    std::shared_ptr<PeriodicMatslise<Scalar>> periodic;
    std::unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction> eigenfunction;

    PeriodicEigenfunctionWrapper(
            std::shared_ptr<PeriodicMatslise<Scalar>> parent,
            std::unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction> eigenfunction) :
            periodic(std::move(parent)), eigenfunction(std::move(eigenfunction)) {
    }

    virtual Eigen::Array<Scalar, 2, 1> operator()(const Scalar &x) const override {
        return (*eigenfunction)(x);
    };

    virtual Eigen::Array<Scalar, Eigen::Dynamic, 2>
    operator()(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const override {
        return (*eigenfunction)(x);
    }
};

template<typename Scalar>
pair<Scalar, Scalar> calculatePeriodicMatchingError(const Y<Scalar, 1, 2> &yError) {
    auto err = yError.y();
    auto dErr = yError.ydE();
    Scalar error = err(0, 0) * err(1, 1) - err(0, 1) * err(1, 0);
    Scalar dError = dErr(0, 0) * err(1, 1) + err(0, 0) * dErr(1, 1) - dErr(0, 1) * err(1, 0) - err(0, 1) * dErr(1, 0);
    return {error, dError};
}

void pyslise(py::module &m) {
    py::class_<AbstractMatslise<double>::Eigenfunction>(m, "Eigenfunction")
            .def("__call__", [](const AbstractMatslise<double>::Eigenfunction &eigenfunction, double x) {
                return eigenfunction(x);
            })
            .def("__call__", [](const AbstractMatslise<double>::Eigenfunction &eigenfunction, const ArrayXd &x)
                    -> Array<double, 2, Dynamic> {
                return eigenfunction(x).transpose();
            });

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
                         -> Array<double, 2, Dynamic> {
                     const Vector2d &right = _right ? *_right : left;
                     Array<double, Dynamic, 2> ys = (*m.eigenfunction(E, make_y(left), make_y(right), index))(xs);
                     return ys.transpose();
                 }, R""""(\
Calculate the eigenfunction corresponding to the eigenvalue E in the points xs.

:param float E: the eigenvalue.
:param (float,float) left, right: the boundary conditions.
:param xs: the points to calculate the eigenfunction for.

:returns: a pair of lists which each a length of len(xs). The first list contains the values of the eigenfunction in the points xs. The second contains the derivative of the eigenfunction in those points.
)"""",
                 py::arg("E"), py::arg("xs"), py::arg("left"), py::arg("right") = optional<Vector2d>(),
                 py::arg("index") = -1)
            .def("eigenfunction",
                 [](AbstractMatslise<double> &m, double E, const Vector2d &left,
                    const optional<Vector2d> &_right, int index)
                         -> std::unique_ptr<AbstractMatslise<double>::Eigenfunction> {
                     const Vector2d &right = _right ? *_right : left;
                     return m.eigenfunction(E, make_y(left), make_y(right), index);
                 }, R""""(\
Calculate the eigenfunction corresponding to the eigenvalue E in the points xs.

:param float E: the eigenvalue.
:param (float,float) left, right: the boundary conditions.
:param xs: the points to calculate the eigenfunction for.

:returns: a pair of lists which each a length of len(xs). The first list contains the values of the eigenfunction in the points xs. The second contains the derivative of the eigenfunction in those points.
)"""",
                 py::arg("E"), py::arg("left"), py::arg("right") = optional<Vector2d>(),
                 py::arg("index") = -1, py::keep_alive<0, 1>())
            .def("eigenfunctionCalculator",
                 [](AbstractMatslise<double> &m, double E, const Vector2d &left,
                    const optional<Vector2d> &_right, int index)
                         -> std::unique_ptr<AbstractMatslise<double>::Eigenfunction> {
                     const Vector2d &right = _right ? *_right : left;
                     return m.eigenfunction(E, make_y(left), make_y(right), index);
                 }, R""""(\
Returns the eigenfunction corresponding to the eigenvalue E as a python function. The returned function can be evaluated in all the points in the domain.

:param float E: the eigenvalue.
:param (float,float) left, [right]: the boundary conditions.
:param int [index]: the index of the eigenvalue

:returns: a function that takes a value and returns a tuple with the eigenfunction and derivative in that value.
)"""",
                 py::arg("E"), py::arg("left"), py::arg("right") = optional<Vector2d>(), py::arg("index") = -1,
                 py::keep_alive<0, 1>())
            .def_property_readonly("domain", [](const AbstractMatslise<double> &matslise) {
                return std::pair<double, double>(matslise.domain.min(), matslise.domain.max());
            });


    py::class_<Matslise<>, AbstractMatslise<double>>(m, "Pyslise", R""""(\
>>> from math import pi, cos
>>> import numpy as np
>>> p = Pyslise(lambda x: 2*cos(2*x), 0, pi, 1e-8)
>>> i, E = p.eigenvaluesByIndex(0, 1, (0,1))[0]
>>> abs(E - -0.1102488054219) < 1e-6
True
>>> f = p.eigenfunction(E, left=(0,1), index=i)
>>> all(abs(f(x)[0] - p.eigenfunction(E, [x], left=(0,1), index=i)[0][0]) < 1e-6 for x in np.linspace(0, pi, 100))
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
                     return make_tuple(y0.y(), theta);
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
                     return make_tuple(y0.y(), y0.ydE(), theta);
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
                return p.sectors[i].get();
            }, py::return_value_policy::reference);

    py::class_<MatsliseHalf<>, AbstractMatslise<double>>(m, "PysliseHalf")
            .def(py::init([](const function<double(double)> &V, double xmax, double tolerance) {
                return new MatsliseHalf<>(V, xmax, tolerance);
            }), R""""(\
In the __init__ function all needed data will be precomputed to effectively solve the Schrödinger equation with given potential on the interval [min; max]. Because of the precomputation the function V is only evaluated at the moment of initialisation. Calling other methods after the object is created never will evaluate V.

Note: only one of steps and tolerance have to be set.

:param (float)->float V: the potential.
:param float min, max: the ends of the domain.
:param int steps: the number of steps to take.
:param int tolerance: automatically choose steps with at least the given accuracy.
)"""", py::arg("V"), py::arg("xmax"), py::arg("tolerance") = 1e-8);

    py::class_<PeriodicMatslise<>, shared_ptr<PeriodicMatslise<>>>(m, "PyslisePeriodic", R""""(\
)"""")
            .def(py::init([](const function<double(double)> &V, double min, double max, double tolerance) {
                return std::make_shared<PeriodicMatslise<>>(V, min, max, tolerance);
            }), R""""(\
:param (float)->float V: the potential.
:param float min, max: the ends of the domain.
:param int tolerance: automatically choose steps with at least the given accuracy.
)"""", py::arg("V"), py::arg("min"), py::arg("max"), py::arg("tolerance") = 1e-8)
            .def("propagate",
                 [](const std::shared_ptr<PeriodicMatslise<>> &m, double E, const Matrix2d &y, double a,
                    double b) -> Matrix2d {
                     return m->propagate(E, make_y(y), a, b).first.y();
                 }, R""""(\
For a given E and initial condition in point a, propagate the solution of the Schrödinger equation to the point b.

:param float E: the fixed eigenvalue to use. This value doesn't have to be a true eigenvalue.
:param (float,float) y: the initial condition in point a.
:param float a: the start of the propagation.
:param float b: the point in which the solution is sought.

:returns: a tuple of length two. The first element is a tuple with the value and the derivative in the point b. The second element is the angle found by the Prüfer transformation.
)"""",
                 py::arg("E"), py::arg("y"), py::arg("a"), py::arg("b"))
            .def("__error",
                 [](const shared_ptr<PeriodicMatslise<>> &m, double E)
                         -> tuple<double, double, Array2d> {
                     auto error = m->matchingError(E);
                     return std::tuple_cat(calculatePeriodicMatchingError<double>(error.first),
                                           std::tuple<Array2d>{error.second});
                 },
                 py::arg("E"))
            .def("eigenpairsByIndex",
                 [](const shared_ptr<PeriodicMatslise<>> &m, int iMin, int iMax)
                         -> vector<tuple<int, double, vector<unique_ptr<AbstractMatslise<double>::Eigenfunction>>>> {
                     vector<tuple<int, double, vector<unique_ptr<AbstractMatslise<double>::Eigenfunction>>>> result;
                     auto pairs = m->eigenpairsByIndex(iMin, iMax);
                     result.reserve(pairs.size());
                     for (auto &eigenpair : pairs) {
                         vector<unique_ptr<AbstractMatslise<double>::Eigenfunction>> eigenfunctions;
                         eigenfunctions.reserve(get<2>(eigenpair).size());
                         for (auto &f : get<2>(eigenpair))
                             eigenfunctions.emplace_back(std::make_unique<PeriodicEigenfunctionWrapper<double>>(
                                     m, std::move(f)));
                         result.emplace_back(get<0>(eigenpair), get<1>(eigenpair), std::move(eigenfunctions));
                     }
                     return result;
                 }, "",
                 py::arg("iMin"), py::arg("iMax"))
            .def("eigenvaluesByIndex", &PeriodicMatslise<>::eigenvaluesByIndex,
                 "", py::arg("iMin"), py::arg("iMax"))
            .def("eigenfunction",
                 [](const shared_ptr<PeriodicMatslise<>> &m, double E)
                         -> vector<unique_ptr<AbstractMatslise<double>::Eigenfunction>> {
                     vector<std::unique_ptr<PeriodicMatslise<>::Eigenfunction>> fs;
                     fs.reserve(2);
                     for (auto &f : m->eigenfunction(E))
                         fs.emplace_back(std::make_unique<PeriodicEigenfunctionWrapper<double>>(m, std::move(f)));
                     return fs;
                 }, "",
                 py::arg("E"));

    py::class_<LiouvilleTransformation<double>>(m, "LiouvilleTransformation")
            .def(py::init([](const function<double(double)> &p, const function<double(double)> &q,
                             const function<double(double)> &w, double min, double max) {
                return LiouvilleTransformation<double>({min, max}, p, q, w);
            }), R""""(\
)"""", py::arg("p"), py::arg("q"), py::arg("w"), py::arg("min"), py::arg("max"))
            .def("V", &LiouvilleTransformation<double>::V, py::arg("x"))
            .def("r2x", &LiouvilleTransformation<double>::r2x, py::arg("r"))
            .def("x2r", &LiouvilleTransformation<double>::x2r, py::arg("x"))
            .def("pwx", &LiouvilleTransformation<double>::pwx, py::arg("x"))
            .def_property_readonly("rDomain", [](const LiouvilleTransformation<double> &lt) {
                return std::pair{lt.domain.min(), lt.domain.max()};
            })
            .def_property_readonly("xDomain", [](const LiouvilleTransformation<double> &lt) {
                auto domain = lt.xDomain();
                return std::pair{domain.min(), domain.max()};
            })
            .def_property_readonly("pieces", [](const LiouvilleTransformation<double> &lt) {
                std::vector<std::pair<double, double>> r;
                r.reserve(lt.pieces.size());
                std::transform(lt.pieces.begin(), lt.pieces.end(), std::back_inserter(r), [](auto p) {
                    return std::pair{p.r.min, p.r.max};
                });
                return r;
            });


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