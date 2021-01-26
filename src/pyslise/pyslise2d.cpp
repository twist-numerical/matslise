#include "module.h"
#include "pybind11/operators.h"

void pyslise2d(py::module &m) {
    py::class_ < Eigenfunction2D < double >> (m, "Eigenfunction2D")
            .def("__call__", [](const Eigenfunction2D<double> &eigenfunction, double x, double y) {
                return eigenfunction(x, y);
            })
            .def("__call__", [](const Eigenfunction2D<double> &eigenfunction, const ArrayXd &x, const ArrayXd &y) {
                return eigenfunction(x, y);
            });

    py::class_ < Eigenfunction2D < double, true >> (m, "Eigenfunction2DWithDerivatives")
            .def("__call__", [](const Eigenfunction2D<double> &eigenfunction, double x, double y) {
                return eigenfunction(x, y);
            })
            .def("__call__", [](const Eigenfunction2D<double> &eigenfunction, const ArrayXd &x, const ArrayXd &y) {
                return eigenfunction(x, y);
            });

    py::class_ < AbstractMatslise2D < double > , shared_ptr < AbstractMatslise2D < double >> > (m, "AbstractPyslise2D")
            .def("eigenvalue", &AbstractMatslise2D<double>::eigenvalue, R""""(\
By using the algorithm of Newton-Raphson the closest eigenvalue around ``start`` will be searched.

:param float start: the initial guess.
:returns: the eigenvalue found starting with ``guess``. Note that the found eigenvalue doesn't necessarily is the closest.
)"""", py::arg("start"))
            .def("eigenvalues", &AbstractMatslise2D<double>::eigenvalues, R""""(\
Calculate the eigenvalues in an interval [Emin; Emax]. The boundary conditions have to be specified.

:param float Emin Emax: the start and end point of the range that will be searched.
:returns: a list of found eigenvalues. When one has a larger multiplicity it is repeated.
)"""", py::arg("Emin"), py::arg("Emax"))
            .def("eigenvaluesByIndex", &AbstractMatslise2D<double>::eigenvaluesByIndex, R""""(\
Calculate all eigenvalues with index between Imin and Imax. The first eigenvalue has index 0. Imin inclusive, Imax exclusive.

:param int Imin: the first eigenvalue to find, by index.
:param int Imax: only the first Imax eigenvalues will be considered.

:returns: a list of eigenvalues.
)"""", py::arg("Imin"), py::arg("Imax"))
            .def("eigenvalueError", &AbstractMatslise2D<double>::eigenvalueError, R""""(\
Estimate an error of a given eigenvalue by using a lower order method.

:param float E: an eigenvalue to estimate an error for.
:returns: the estimated error for the given eigenvalue.
)"""", py::arg("E"))
            .def("eigenfunction",
                 [](const Matslise2D<> &se2d, double E) -> vector <Eigenfunction2D<double>> {
                     vector <Eigenfunction2D<double>> result;
                     for (auto &f : se2d.eigenfunction(E))
                         result.push_back(f);
                     return result;
                 }, R""""(\
Returns a list if eigenfunctions corresponding to the eigenvalue E as python functions. The returned functions can be evaluated in all the points in the domain.

:param float E: the eigenvalue.

:returns: a list of functions (depending on multiplicity) each taking a x-value and a y-value (or a two lists to evaluate a grid) and returning the value of that eigenfunction in (x, y).
)"""", py::arg("E"))
            .def("eigenfunctionWithDerivatives",
                 [](const Matslise2D<> &se2d, double E) -> vector <Eigenfunction2D<double, true>> {
                     vector <Eigenfunction2D<double, true>> result;
                     for (auto &f : se2d.eigenfunctionWithDerivatives(E))
                         result.push_back(f);
                     return result;
                 }, R""""(\
Returns a list if eigenfunctions corresponding to the eigenvalue E as python functions. The returned functions can be evaluated in all the points in the domain.

:param float E: the eigenvalue.

:returns: a list of functions (depending on multiplicity) each taking a x-value and a y-value (or a two lists to evaluate a grid) and returning the value, the x-derivative and the y-derivative of that eigenfunction in (x, y).
)"""", py::arg("E"))
            .def("estimateIndex", &AbstractMatslise2D<double>::estimateIndex, py::arg("E"));


    py::class_ < Matslise2D<>, AbstractMatslise2D < double >, shared_ptr < Matslise2D < double >> > (m, "Pyslise2D", R""""(\
>>> p_ixaru = Pyslise2D(lambda x, y: (1+x*x)*(1+y*y), -5.5,5.5, -5.5,5.5)
>>> [(int(i), round(E, 4), int(m)) for i, E, m in p_ixaru.eigenvaluesByIndex(0, 4)]
[(0, 3.1959, 1), (1, 5.5267, 2), (3, 7.5578, 1), (4, 8.0313, 1)]
>>> p_henon = Pyslise2D(lambda x, y: (x*x + y*y) + 2*0.11180339887*y*(x*x-y*y/3), -11,11, -11,11, x_symmetric=True)
>>> [(int(i), round(E, 4), int(m)) for i, E, m in p_henon.eigenvaluesByIndex(0, 4)]
[(0, 1.9972, 1), (1, 3.9802, 2), (3, 5.9125, 1), (4, 5.9707, 2)]
)"""")
            .def(py::init([](const function<double(double, double)> &V,
                             double xmin, double xmax, double ymin, double ymax,
                             bool symmetric,
                             int x_count, double x_tol,
                             int y_count, double y_tol,
                             double tol, int N, int stepsPerSector) {
                     py::gil_scoped_release release;
                     if (x_count != -1 && x_tol != -1) {
                         throw invalid_argument("Not both 'x_count' and 'x_tolerance' can be set.");
                     }
                     if (x_count == -1 && x_tol == -1) {
                         if (tol != -1)
                             x_tol = tol;
                         else
                             throw invalid_argument("One of 'tolerance', 'x_count' or 'x_tolerance' must be set.");
                     }
                     if (y_count != -1 && y_tol != -1) {
                         throw invalid_argument("Not both 'y_count' and 'y_tolerance' can be set.");
                     }
                     if (y_count == -1 && y_tol == -1) {
                         if (tol != -1)
                             y_tol = tol;
                         else
                             throw invalid_argument("One of 'tolerance', 'y_count' or 'y_tolerance' must be set.");
                     }
                     Matslise2D<>::Config config;
                     config.tolerance = tol;
                     config.basisSize = N;
                     config.xSymmetric = symmetric;
                     config.stepsPerSector = stepsPerSector;

                     if (x_count != -1)
                         config.xSectorBuilder = sector_builder::uniform<Matslise<>>(x_count);
                     else
                         config.xSectorBuilder = sector_builder::automatic<Matslise<>>(x_tol);

                     if (y_count != -1)
                         config.ySectorBuilder = sector_builder::uniform<Matslise2D<>>(y_count);
                     else
                         config.ySectorBuilder = sector_builder::automatic<Matslise2D<>>(y_tol);

                     return make_unique<Matslise2D<>>([V](double x, double y) -> double {
                         py::gil_scoped_acquire acquire;
                         return V(x, y);
                     }, Rectangle < double, 2 > {xmin, xmax, ymin, ymax}, config);
                 }),
                 R""""(\
In the __init__ function all needed data will be precomputed to effectively solve the given Schrödinger equation on the domain. Because of the precomputation the function V is only evaluated at the moment of initalisation. Calling other methods when the object is created will never evaluate V.

Note: steps along the y-axis are more computational expensive.

:param (float,float)->float V: the potential of the Schrödinger equation to solve
:param float xmin, xmax, ymin, ymax: the domain to work on.
:param bool x_symmetric: If this is true the potential is assumed to be even in the x-direction.
:param float x_count, x_tolerance: Use only one of these. This is a guidance on how to pick the number of steps along the x-axis. With x_count pyslise will make uniform_steps. With x_tolerance the steps will be picked to try to keep the error lower than x_tolerance.
:param float y_count, y_tolerance: analogous the x_count and y_tolerance.
:param float tolerance: if none of x_count, x_tolerance, y_count or y_tolerance. x_tolerance and y_tolerance will be set to tolerance. Defaults to 1e-7.

The next set of parameters are more advanced. Tweaking these can be useful when the required accuracy isn't reached.

:param int N: the number of used basis functions on each sector. Defaults to 12.
:param int steps_per_sector: the number of steps that will be taken per sector (in de y-direction). Defaults to 2.
)"""",
                 py::arg("V"),
                 py::arg("xmin"), py::arg("xmax"), py::arg("ymin"), py::arg("ymax"),
                 py::arg("x_symmetric") = false,
                 py::arg("x_count") = -1, py::arg("x_tolerance") = -1,
                 py::arg("y_count") = -1, py::arg("y_tolerance") = -1,
                 py::arg("tolerance") = 1e-7, py::arg("N") = 12, py::arg("steps_per_sector") = 2)
            .def("__matchingError", [](const Matslise2D<> &se2d, double const &E) -> pair<double, double> {
                return se2d.matchingError(E);
            }, R""""(\
Compute the error given a guess for E. This error is the result of the requirement that the found eigenfunctions are continues. The error expresses how 'discontinues' the corresponding eigenfunction would be.

:param float E: the guessed eigenvalue.
:returns: A tuple with the computed error and the derivative of that error with respect to E.
)"""", py::arg("E"))
            .def("__matchingErrors", [](const Matslise2D<> &se2d, double E) -> vector <pair<double, double>> {
                return se2d.matchingErrors(E);
            }, R""""(\
Just like Pyslise2D::matchingError(E) computes this function the discontinuity of the eigenfunction. The corresponding eigenfunction will be continuous once any of the N returned values is zero.

:param float E: the guessed eigenvalue.
:returns: A list of tuples with each of the computed errors and its derivative with respect to E.
)"""", py::arg("E"))
            .def("__propagate",
                 [](const Matslise2D<> &m, double E, const MatrixXd &y, const MatrixXd &dy, double a, double b) ->
                         pair <MatrixXd, MatrixXd> {
                     Y<double, Dynamic> y0(m.config.basisSize);
                     y0.block() = y;
                     y0.block(dX) = dy;
                     return unpackY(m.propagate(E, y0, a, b)).first;
                 },
                 py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
            .def_property_readonly("__N", [](const Matslise2D<> &m) {
                return m.config.basisSize;
            }, "The number of basis functions used on each sector")
            .def_property_readonly("__M", [](Matslise2D<> &p) -> vector <MatrixXd> * {
                auto l = new vector<MatrixXd>(static_cast<vector<MatrixXd>::size_type>(p.sectors.size() - 1));
                for (unsigned long i = 0; i < p.sectors.size() - 1; ++i)
                    l->at(i) = p.M[i];
                return l;
            })
            .def_property_readonly("__matchpoint", [](const Matslise2D<> &p) -> double {
                return p.sectors[p.matchIndex]->max;
            })
            .def_property_readonly("__sectors", [](const Matslise2D<> &p) {
                return p.sectors;
            });

    py::class_<Matslise2DHalf<>, AbstractMatslise2D<double>, shared_ptr<Matslise2DHalf<double>>>(m, "Pyslise2DHalf")
            .def(py::init([](const function<double(double, double)> &V,
                             double xmin, double xmax, double ymax,
                             bool symmetric,
                             int x_count, double x_tol,
                             int y_count, double y_tol,
                             double tol, int N, int stepsPerSector) {
                     py::gil_scoped_release release;
                     if (x_count != -1 && x_tol != -1) {
                         throw invalid_argument("Not both 'x_count' and 'x_tolerance' can be set.");
                     }
                     if (x_count == -1 && x_tol == -1) {
                         if (tol != -1)
                             x_tol = tol;
                         else
                             throw invalid_argument("One of 'tolerance', 'x_count' or 'x_tolerance' must be set.");
                     }
                     if (y_count != -1 && y_tol != -1) {
                         throw invalid_argument("Not both 'y_count' and 'y_tolerance' can be set.");
                     }
                     if (y_count == -1 && y_tol == -1) {
                         if (tol != -1)
                             y_tol = tol;
                         else
                             throw invalid_argument("One of 'tolerance', 'y_count' or 'y_tolerance' must be set.");
                     }
                     Matslise2D<>::Config config;
                     config.tolerance = tol;
                     config.basisSize = N;
                     config.xSymmetric = symmetric;
                     config.stepsPerSector = stepsPerSector;

                     if (x_count != -1)
                         config.xSectorBuilder = sector_builder::uniform<Matslise<>>(x_count);
                     else
                         config.xSectorBuilder = sector_builder::automatic<Matslise<>>(x_tol);

                     if (y_count != -1)
                         config.ySectorBuilder = sector_builder::uniform<Matslise2D<>>(y_count);
                     else
                         config.ySectorBuilder = sector_builder::automatic<Matslise2D<>>(y_tol);

                     return make_unique<Matslise2DHalf<>>([V](double x, double y) -> double {
                         py::gil_scoped_acquire acquire;
                         return V(x, y);
                     }, Rectangle<double, 2>{xmin, xmax, -ymax, ymax}, config);
                 }),
                 R""""(\
In the __init__ function all needed data will be precomputed to effectively solve the given Schrödinger equation on the domain. Because of the precomputation the function V is only evaluated at the moment of initalisation. Calling other methods when the object is created will never evaluate V.

Note: steps along the y-axis are more computational expensive.

:param (float,float)->float V: the potential of the Schrödinger equation to solve
:param float xmin, xmax, ymax: the domain to work on.
:param bool x_symmetric: If this is true the potential is assumed to be even in the x-direction.
:param float x_count, x_tolerance: Use only one of these. This is a guidance on how to pick the number of steps along the x-axis. With x_count pyslise will make uniform_steps. With x_tolerance the steps will be picked to try to keep the error lower than x_tolerance.
:param float y_count, y_tolerance: analogous the x_count and y_tolerance.
:param float tolerance: if none of x_count, x_tolerance, y_count or y_tolerance. x_tolerance and y_tolerance will be set to tolerance. Defaults to 1e-7.

The next set of parameters are more advanced. Tweaking these can be useful when the required accuracy isn't reached.

:param int N: the number of used basis functions on each sector. Defaults to 12.
:param int steps_per_sector: the number of steps that will be taken per sector (in de y-direction). Defaults to 2.
)"""",
                         py::arg("V"),
                         py::arg("xmin"), py::arg("xmax"),
                         py::arg("ymax"),
                         py::arg("x_symmetric") = false,
                         py::arg("x_count") = -1,
                         py::arg("x_tolerance") = -1,
                         py::arg("y_count") = -1,
                         py::arg("y_tolerance") = -1,
                         py::arg("tolerance") = 1e-7,
                         py::arg("N") = 12,
                         py::arg("steps_per_sector") = 2);

    py::class_<Matslise2D<>::Sector, std::unique_ptr<Matslise2D<>::Sector, py::nodelete>>(m, "Pyslise2DSector")
            .def_property_readonly("eigenvalues", [](Matslise2D<>::Sector &s) -> vector<double> * {
                auto l = new vector<double>(s.se2d->config.basisSize);
                for (unsigned long i = 0; i < l->size(); ++i)
                    l->at(i) = s.eigenvalues[i];
                return l;
            })
            .def("error", &Matslise2D<>::Sector::error)
            .def("estimateIndex",
                 [](const Matslise2D<>::Sector &sector, double E, const MatrixXd &y, const MatrixXd &dy) -> Index {
                     Y<double, Dynamic> y0(sector.se2d->config.basisSize);
                     y0.block() = y;
                     y0.block(dX) = dy;
                     return sector.propagateWithIndex(E, y0).second;
                 },
                 py::arg("E"), py::arg("y"), py::arg("dy"))
            .def_property_readonly("matslise",
                                   [](Matslise2D<>::Sector &s) -> AbstractMatslise<double> * {
                                       return s.matslise.get();
                                   }, py::return_value_policy::reference_internal)
            .def_readonly("matscs", &Matslise2D<>::Sector::matscs, py::return_value_policy::reference)
            .def_readonly("min", &Matslise2D<>::Sector::min)
            .def_readonly("max", &Matslise2D<>::Sector::max);
}