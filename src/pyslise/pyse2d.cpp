#include "module.h"

void pySE2d(py::module &m) {
    py::class_<SE2D<>>(m, "PySE2d")
            .def(py::init([](function<double(double, double)> V,
                             double xmin, double xmax, double ymin, double ymax,
                             bool symmetric,
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
                             .nested(o1.symmetric(symmetric));
                     if (y_count != -1)
                         o2.sectorCount(y_count);
                     else
                         o2.tolerance(y_tol);
                     return unique_ptr<SE2D<>>(new SE2D<>(V, {{xmin, xmax}, ymin, ymax}, o2));
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
                 py::arg("symmetric") = true,
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
            .def("eigenfunction", &SE2D<>::eigenfunctionCalculator, R""""(\
Returns a list if eigenfunctions corresponding to the eigenvalue E as python functions. The returned functions can be evaluated in all the points in the domain.

:param float E: the eigenvalue.

:returns: a list of function that takes a x-value and a y-value and returns the value of the eigenfunction in (x, y).
)"""",
                 py::arg("E"))
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
            .def("firstEigenvalue", &SE2D<>::findFirstEigenvalue)
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
            .def_property_readonly("__sectors", [](SE2D<> &p) -> std::vector<SE2D<>::Sector *> {
                vector<SE2D<>::Sector *> l(p.sectorCount);
                for (int i = 0; i < p.sectorCount; ++i)
                    l[i] = p.sectors[i];
                return l;
            });

    py::class_<SE2D<>::Sector, std::unique_ptr<SE2D<>::Sector, py::nodelete>>(m, "PySE2dSector")
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
}