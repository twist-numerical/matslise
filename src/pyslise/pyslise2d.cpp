#include "module.h"

void pyslise2d(py::module &m) {
    py::class_<AbstractMatslise2D<double>, shared_ptr<AbstractMatslise2D<double>>>(m, "AbstractPyslise2D")
            .def("eigenvalue", &AbstractMatslise2D<double>::eigenvalue, R""""(\
By using the algorithm of Newton-Raphson the closest eigenvalue around ``start`` will be searched. It keeps executing this algorithm until either the number of iterations is reached or the error drops below tolerance.

:param float start: the initial guess.
:param float tolerance: one of the stopping conditions. Defaults to 1e-9.
:param int iterations: the other stopping conditions. Defaults to 30.
:param float min_tolerance: if the maximum number of iterations is reached and the error is smaller than ``min_tolerance`` then the found value will be counted is eigenvalue. Defaults to 1e-5.
:returns: the eigenvalue found starting with ``guess``. Note that the found eigenvalue doesn't necessarily is the closest .
)"""", py::arg("start"))
            .def("eigenvalues", &AbstractMatslise2D<double>::eigenvalues, R""""(\
This heuristics tries to find all the eigenvalues within a certain range [Emin, Emax]. Because this heuristics isn't an algorithm, it is certainly not certain that all eigenvalues are found. In short: the heuristics starts with a few initial guesses and tries to find all eigenvalues that it can 'see' from that first guess.

It is not a good idea to make the number of initial values large. This will increase computation time and, more importantly, it won't be necessarily better.

:param float Emin Emax: the start and end point of the range that will be searched.
:param int initial_values: the number of starting guesses that will be used. Defaults to 16.
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
                 [](const AbstractMatslise2D<double> &se2d, double E, const ArrayXd &x, const ArrayXd &y)
                         -> vector<ArrayXXd> {
                     vector<ArrayXXd> result;
                     for (auto &f : se2d.eigenfunction(E))
                         result.push_back(f(x, y));
                     return result;
                 }, R""""(\
Compute all the corresponding eigenfunctions for a given eigenvalue. Most of the time this will return a singleton list. But it is possible that this eigenvalue has a higher multiplicity, so more eigenfunctions will be returned. On the other hand, when the given value for E isn't an eigenvalue then there doesn't exist an eigenfunction, so the returned list will be empty.

:param float E: the eigenvalue to compute eigenfunctions for.
:param [float] x y: the x and y values of the points to evaluate the eigenfunctions in.
:returns: a list of len(x) by len(y) grids of values. In each grid the value on position i, j is that eigenfunction evaluated in point x[i], y[j].
)"""", py::arg("E"), py::arg("x"), py::arg("y"))
            .def("eigenfunction",
                 [](const Matslise2D<> &se2d, double E) -> vector<function<double(double, double)>> {
                     vector<function<double(double, double)>> result;
                     for (auto &f : se2d.eigenfunction(E))
                         result.push_back(f);
                     return result;
                 }, R""""(\
Returns a list if eigenfunctions corresponding to the eigenvalue E as python functions. The returned functions can be evaluated in all the points in the domain.

:param float E: the eigenvalue.

:returns: a list of functions (depending on multiplicity) each taking a x-value and a y-value and returning the value of that eigenfunction in (x, y).
)"""", py::arg("E"))
            .def("eigenfunctionDerivatives",
                 [](const AbstractMatslise2D<double> &se2d, double E, const ArrayXd &x, const ArrayXd &y)
                         -> vector<tuple<ArrayXXd, ArrayXXd, ArrayXXd>> {
                     vector<tuple<ArrayXXd, ArrayXXd, ArrayXXd>> result;
                     for (auto &f : se2d.eigenfunctionWithDerivatives(E))
                         result.push_back(f(x, y));
                     return result;
                 }, R""""(\
Compute all the corresponding eigenfunctions for a given eigenvalue. Most of the time this will return a singleton list. But it is possible that this eigenvalue has a higher multiplicity, so more eigenfunctions will be returned. On the other hand, when the given value for E isn't an eigenvalue then there doesn't exist an eigenfunction, so the returned list will be empty.

:param float E: the eigenvalue to compute eigenfunctions for.
:param [float] x y: the x and y values of the points to evaluate the eigenfunctions in.
:returns: a list of three len(x) by len(y) grids of values. For each triplet of grids, position i, j contains the value (resp. x-derivative and y-derivative) of that eigenfunction evaluated in the point x[i], y[j].
)"""", py::arg("E"), py::arg("x"), py::arg("y"))
            .def("eigenfunctionDerivatives",
                 [](const Matslise2D<> &se2d, double E) -> vector<function<
                         tuple<double, double, double>(double, double)>> {
                     vector<function<
                             tuple<double, double, double>(double, double)>> result;
                     for (auto &f : se2d.eigenfunctionWithDerivatives(E))
                         result.push_back(f);
                     return result;
                 }, R""""(\
Returns a list if eigenfunctions corresponding to the eigenvalue E as python functions. The returned functions can be evaluated in all the points in the domain.

:param float E: the eigenvalue.

:returns: a list of functions (depending on multiplicity) each taking a x-value and a y-value and returning the value, the x-derivative and the y-derivative of that eigenfunction in (x, y).
)"""", py::arg("E"))
            .def("estimateIndex", &AbstractMatslise2D<double>::estimateIndex, py::arg("E"));


    py::class_<Matslise2D<>, AbstractMatslise2D<double>, shared_ptr<Matslise2D<double>>>(m, "Pyslise2D")
            .def(py::init([](const function<double(double, double)> &V,
                             double xmin, double xmax, double ymin, double ymax,
                             bool symmetric,
                             int x_count, double x_tol,
                             int y_count, double y_tol,
                             double tol, int N, int stepsPerSector) {
                     py::gil_scoped_release release;
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
                     }, Rectangle<double, 2>{xmin, xmax, ymin, ymax}, config);
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
                 py::arg("symmetric") = false,
                 py::arg("x_count") = -1, py::arg("x_tolerance") = -1,
                 py::arg("y_count") = -1, py::arg("y_tolerance") = -1,
                 py::arg("tolerance") = -1, py::arg("N") = 12, py::arg("steps_per_sector") = 3)
            .def("matchingError", [](const Matslise2D<> &se2d, double const &E) -> pair<double, double> {
                return se2d.matchingError(E);
            }, R""""(\
Compute the error given a guess for E. This error is the result of the requirement that the found eigenfunctions are continues. The error expresses how 'discontinues' the corresponding eigenfunction would be.

:param float E: the guessed eigenvalue.
:returns: A tuple with the computed error and the derivative of that error with respect to E.
)"""", py::arg("E"))
            .def("matchingErrors", [](const Matslise2D<> &se2d, double E) -> vector<pair<double, double>> {
                return se2d.matchingErrors(E);
            }, R""""(\
Just like Pyslise2D::calculateError(E) computes this function the discontinuity of the eigenfunction. The corresponding eigenfunction will be continuous once any of the N returned values is zero.

:param float E: the guessed eigenvalue.
:returns: A list of tuples with each of the computed errors and its derivative with respect to E.
)"""", py::arg("E"))
            .def("__propagate",
                 [](const Matslise2D<> &m, double E, const MatrixXd &y, const MatrixXd &dy, double a, double b) ->
                         pair<MatrixXd, MatrixXd> {
                     Y<double, Dynamic> y0(m.config.basisSize);
                     y0.block() = y;
                     y0.block(dX) = dy;
                     return unpackY(m.propagate(E, y0, a, b)).first;
                 },
                 py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
            .def_property_readonly("__N", [](const Matslise2D<> &m) {
                return m.config.basisSize;
            }, "The number of basis functions used on each sector")
            .def_property_readonly("__M", [](Matslise2D<> &p) -> vector<MatrixXd> * {
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
:param float x_count, x_tolerance: Use only one of these. This is a guidance on how to pick the number of steps along the x-axis. With x_count pyslise will make uniform_steps. With x_tolerance the steps will be picked to try to keep the error lower than x_tolerance.
:param float y_count, y_tolerance: analogous the x_count and y_tolerance.
:param float tolerance: if none of x_count, x_tolerance, y_count or y_tolerance. x_tolerance and y_tolerance will be set to tolerance

The next set of parameters are more advanced and can be useful to tweak when the required accuracy isn't reached.

:param int N: the number of used basis functions on each sector. Defaults to 12.
:param int in_sector_count: the number of steps that will be taken per sector (in de y-direction). Defaults to 2.
:param in grid_points: the number of points that will be used to calculate the quadratures. Defaults to 60.
)"""",
                 py::arg("V"),
                 py::arg("xmin"), py::arg("xmax"), py::arg("ymax"),
                 py::arg("symmetric") = false,
                 py::arg("x_count") = -1, py::arg("x_tolerance") = -1,
                 py::arg("y_count") = -1, py::arg("y_tolerance") = -1,
                 py::arg("tolerance") = -1, py::arg("N") = 12, py::arg("steps_per_sector") = 3);

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