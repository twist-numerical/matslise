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
)"""", py::arg("E"))
            .def("eigenvalues", &AbstractMatslise3D<double>::eigenvalues, R""""(\
This heuristics tries to find all the eigenvalues within a certain range [Emin, Emax]. Because this heuristics isn't an algorithm, it is certainly not certain that all eigenvalues are found. In short: the heuristics starts with a few initial guesses and tries to find all eigenvalues that it can 'see' from that first guess.

It is not a good idea to make the number of initial values large. This will increase computation time and, more importantly, it won't be necessarily better.

:param float Emin Emax: the start and end point of the range that will be searched.
:param int initial_values: the number of starting guesses that will be used. Defaults to 16.
:returns: a list of found eigenvalues. When one has a larger multiplicity it is repeated.
)"""", py::arg("Emin"), py::arg("Emax"))
            .def("eigenvaluesByIndex", &AbstractMatslise3D<double>::eigenvaluesByIndex, R""""(\
Calculate all eigenvalues with index between Imin and Imax. The first eigenvalue has index 0. Imin inclusive, Imax exclusive.

:param int Imin: the first eigenvalue to find, by index.
:param int Imax: only the first Imax eigenvalues will be considered.

:returns: a list of eigenvalues.
)"""", py::arg("Imin"), py::arg("Imax"))
            .def("estimateIndex", &AbstractMatslise3D<double>::estimateIndex, py::arg("E"));

    py::class_<Matslise3D<double>, AbstractMatslise3D<double>, shared_ptr<Matslise3D<double>>>(m, "Pyslise3D")
            .def(py::init([](const function<double(double, double, double)> &V,
                             double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                             double tol,
                             int x_count, double x_tol,
                             int y_count, double y_tol,
                             int z_count, double z_tol,
                             bool x_symmetric, bool y_symmetric,
                             int x_basis_size, int xy_basis_size,
                             int y_steps_per_sector, int z_steps_per_sector) {
                     py::gil_scoped_release release;
                     if (x_count != -1 && x_tol != -1) {
                         throw invalid_argument("Not both 'x_count' and 'x_tolerance' can be set.");
                     }
                     if (x_count == -1 && x_tol == -1) {
                         if (tol != -1)
                             x_tol = tol;
                         else
                             throw invalid_argument("One of 'x_count' and 'x_tolerance' must be set.");
                     }
                     if (y_count != -1 && y_tol != -1) {
                         throw invalid_argument("Not both 'y_count' and 'y_tolerance' can be set.");
                     }
                     if (y_count == -1 && y_tol == -1) {
                         if (tol != -1)
                             y_tol = tol;
                         else
                             throw invalid_argument("One of 'y_count' and 'y_tolerance' must be set.");
                     }
                     if (z_count != -1 && z_tol != -1) {
                         throw invalid_argument("Not both 'z_count' and 'z_tolerance' can be set.");
                     }
                     if (z_count == -1 && z_tol == -1) {
                         if (tol != -1)
                             z_tol = tol;
                         else
                             throw invalid_argument("One of 'z_count' and 'z_tolerance' must be set.");
                     }
                     Matslise3D<>::Config config;
                     config.tolerance = tol;
                     config.xBasisSize = x_basis_size;
                     config.xyBasisSize = xy_basis_size;
                     config.xSymmetric = x_symmetric;
                     config.ySymmetric = y_symmetric;
                     config.yStepsPerSector = y_steps_per_sector;
                     config.zStepsPerSector = z_steps_per_sector;

                     if (x_count != -1)
                         config.xSectorBuilder = sector_builder::uniform<Matslise<>>(x_count);
                     else
                         config.xSectorBuilder = sector_builder::automatic<Matslise<>>(x_tol);

                     if (y_count != -1)
                         config.ySectorBuilder = sector_builder::uniform<Matslise2D<>>(y_count);
                     else
                         config.ySectorBuilder = sector_builder::automatic<Matslise2D<>>(y_tol);

                     if (z_count != -1)
                         config.zSectorBuilder = sector_builder::uniform<Matslise3D<>>(z_count);
                     else
                         config.zSectorBuilder = sector_builder::automatic<Matslise3D<>>(z_tol);

                     return make_unique<Matslise3D<>>([V](double x, double y, double z) -> double {
                         py::gil_scoped_acquire acquire;
                         return V(x, y, z);
                     }, Rectangle<double, 3>{xmin, xmax, ymin, ymax, zmin, zmax}, config);
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
                 py::arg("x_count") = -1, py::arg("x_tolerance") = -1,
                 py::arg("y_count") = -1, py::arg("y_tolerance") = -1,
                 py::arg("z_count") = -1, py::arg("z_tolerance") = -1,
                 py::arg("x_symmetric") = false, py::arg("y_symmetric") = false,
                 py::arg("x_basis_size") = 12, py::arg("xy_basis_size") = 12,
                 py::arg("y_steps_per_sector") = 3, py::arg("z_steps_per_sector") = 3)
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
            .def_property_readonly("__sectors", [](const Matslise3D<> &p) -> vector<Matslise3D<>::Sector *> {
                vector<Matslise3D<>::Sector *> v;
                for (auto &s : p.sectors)
                    v.emplace_back(s.get());
                return v;
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
