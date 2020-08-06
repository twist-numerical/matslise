#include "module.h"

void pyslise3d(py::module &m) {
    py::class_<AbstractMatslise3D<double>, shared_ptr<AbstractMatslise3D<double>>>(m, "AbstractPyslise3D")
            .def("eigenvalue", &AbstractMatslise3D<double>::eigenvalue, R""""(\
By using the algorithm of Newton-Raphson the closest eigenvalue around ``start`` will be searched. It keeps executing this algorithm until either the number of iterations is reached or the error drops below tolerance.

:param float start: the initial guess.
:param float tolerance: one of the stopping conditions. Defaults to 1e-9.
:param int iterations: the other stopping conditions. Defaults to 30.
:param float min_tolerance: if the maximum number of iterations is reached and the error is smaller than ``min_tolerance`` then the found value will be counted is eigenvalue. Defaults to 1e-5.
:returns: the eigenvalue found starting with ``guess``. Note that the found eigenvalue doesn't necessarily is the closest .
)"""", py::arg("start"))
            .def("eigenvalueError", &AbstractMatslise3D<double>::eigenvalueError, R""""(\
Estimate an error of a given eigenvalue by using a lower order method.

:param float E: an eigenvalue to estimate an error for.
:returns: the estimated error for the given eigenvalue.
)"""", py::arg("E"));

    py::class_<Matslise3D<double>, AbstractMatslise3D<double>, shared_ptr<Matslise3D<double>>>(m, "Pyslise3D")
            .def(py::init([](const function<double(double, double, double)> &V,
                             double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                             double tolerance, int z_count, double z_tol) {
                     py::gil_scoped_release release;
                     if (z_count != -1 && z_tol != -1) {
                         throw invalid_argument("Not both 'y_count' and 'y_tol' can be set.");
                     }
                     if (z_count == -1 && z_tol == -1) {
                         z_tol = tolerance;
                     }
                     SectorBuilder<Matslise3D<double>> sectorBuilder =
                             z_count == -1 ? sector_builder::automatic<Matslise3D<double>>(z_tol)
                                           : sector_builder::uniform<Matslise3D<double>>(z_count);
                     return make_unique<Matslise3D<double>>([V](double x, double y, double z) -> double {
                         py::gil_scoped_acquire acquire;
                         return V(x, y, z);
                     }, Rectangle<3, double>{{{xmin, xmax}, ymin, ymax}, zmin, zmax}, sectorBuilder, tolerance);
                 }),
                 R""""(\
In the __init__ function all needed data will be precomputed to effectively solve the given Schrödinger equation on the domain. Because of the precomputation the function V is only evaluated at the moment of initalisation. Calling other methods when the object is created will never evaluate V.

Note: steps along the y-axis are more computational expensive.

TODO!

:param (float,float)->float V: the potential of the Schrödinger equation to solve
:param float xmin, xmax, ymin, ymax, zmin, zmax: the domain to work on.
:param float tolerance: if none of x_count, x_tolerance, y_count or y_tolerance. x_tolerance and y_tolerance will be set to tolerance
:param float z_count, z_tolerance: Use only one of these. This is a guidance on how to pick the number of steps along the x-axis. With x_count pyslise will make uniform_steps. With x_tolerance the steps will be picked to try to keep the error lower than x_tolerance.
)"""",
                 py::arg("V"),
                 py::arg("xmin"), py::arg("xmax"), py::arg("ymin"), py::arg("ymax"), py::arg("zmin"), py::arg("zmax"),
                 py::arg("tolerance"),
                 py::arg("z_count") = -1, py::arg("z_tolerance") = -1)
            .def("matchingError", [](const Matslise3D<> &s, double const &E) -> pair<double, double> {
                return s.matchingError(E);
            }, R""""(\
Compute the error given a guess for E. This error is the result of the requirement that the found eigenfunctions are continues. The error expresses how 'discontinues' the corresponding eigenfunction would be.

:param float E: the guessed eigenvalue.
:returns: A tuple with the computed error and the derivative of that error with respect to E.
)"""", py::arg("E"))
            .def("matchingErrors", [](const Matslise3D<> &s, double E) -> vector<pair<double, double>> {
                return s.matchingErrors(E);
            }, R""""(\
Just like Pyslise3D::calculateError(E) computes this function the discontinuity of the eigenfunction. The corresponding eigenfunction will be continuous once any of the N returned values is zero.

:param float E: the guessed eigenvalue.
:returns: A list of tuples with each of the computed errors and its derivative with respect to E.
)"""", py::arg("E"))
            .def_property_readonly("__M", [](Matslise3D<> &p) -> vector<MatrixXd> {
                return p.M;
            })
            .def_property_readonly("__sectors", [](Matslise3D<> &p) -> std::vector<Matslise3D<>::Sector *> {
                return p.sectors;
            });

    py::class_<Matslise3D<>::Sector, std::unique_ptr<Matslise3D<>::Sector, py::nodelete>>(m, "Pyslise3DSector")
            .def_property_readonly("eigenvalues", [](Matslise3D<>::Sector &s) -> vector<double> {
                return s.eigenvalues;
            })
            .def_property_readonly(
                    "matslise2d", [](const Matslise3D<>::Sector &s) -> shared_ptr<AbstractMatslise2D<double>> {
                        return s.matslise2d;
                    })
            .def_property_readonly("eigenfunctions", [](const Matslise3D<>::Sector &s) -> vector<ArrayXXd> {
                return s.eigenfunctions_grid;
            })
            .def_readonly("matscs", &Matslise3D<>::Sector::matscs, py::return_value_policy::reference)
            .def_readonly("min", &Matslise3D<>::Sector::min)
            .def_readonly("max", &Matslise3D<>::Sector::max);
}
