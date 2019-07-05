#include "module.h"

void pySlise(py::module &m) {

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
            .def("propagate",
                 [](Matslise<> &m, double E, const Vector2d &y, double a, double b) ->
                         tuple<Vector2d, double> {
                     Y<> y0;
                     double theta;
                     tie(y0, theta) = m.propagate(E, make_y(y), a, b);
                     return make_tuple(y0.y, theta);
                 },
                 "Propagate without dy",
                 py::arg("E"), py::arg("y"), py::arg("a"), py::arg("b"))
            .def("propagate",
                 [](Matslise<> &m, double E, const Vector2d &y, const Vector2d &dy, double a, double b)
                         -> tuple<Vector2d, Vector2d, double> {
                     Y<> y0;
                     double theta;
                     tie(y0, theta) = m.propagate(E, Y<>(y, dy), a, b);
                     return make_tuple(y0.y, y0.dy, theta);
                 },
                 "Propagate without dy",
                 py::arg("E"), py::arg("y"), py::arg("dy"), py::arg("a"), py::arg("b"))
            .def("eigenvalues",
                 [](Matslise<> &m, double Emin, double Emax, const Vector2d &left, const Vector2d &right)
                         -> vector<pair<int, double>> * {
                     return m.computeEigenvalues(Emin, Emax, make_y(left), make_y(right));
                 },
                 py::arg("Emin"), py::arg("Emax"), py::arg("left"), py::arg("right"))
            .def("eigenvaluesByIndex",
                 [](Matslise<> &m, int Imin, int Imax, const Vector2d &left, const Vector2d &right)
                         -> vector<pair<int, double>> * {
                     return m.computeEigenvaluesByIndex(Imin, Imax, make_y(left), make_y(right));
                 },
                 py::arg("Imin"), py::arg("Imax"), py::arg("left"), py::arg("right"))
            .def("eigenvalueError",
                 [](Matslise<> &m, double E, const Vector2d &left, const Vector2d &right)
                         -> double {
                     return m.computeEigenvalueError(E, make_y(left), make_y(right));
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"))
            .def("eigenfunction",
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
            .def("eigenfunction",
                 [](Matslise<> &m, double E, const Vector2d &left,
                    const Vector2d &right) -> std::function<pair<double, double>(double)> {
                     std::function<Y<>(double)> calculator = m.eigenfunctionCalculator(E, make_y(left), make_y(right));
                     return [calculator](double x) -> pair<double, double> {
                         Y<> y = calculator(x);
                         return make_pair(y.y[0], y.y[1]);
                     };
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"))
            .def("error",
                 [](Matslise<> &m, double E, const Vector2d &left, const Vector2d &right)
                         -> tuple<double, double, double> {
                     return m.calculateError(E, make_y(left), make_y(right));
                 },
                 py::arg("E"), py::arg("left"), py::arg("right"))
            .def_readonly("__sectorCount", &Matslise<>::sectorCount)
            .def_readonly("__match", &Matslise<>::match)
            .def_readonly("__min", &Matslise<>::xmin)
            .def_readonly("__max", &Matslise<>::xmax)
            .def("__sector", [](Matslise<> &p, int i) -> Matslise<>::Sector * {
                return p.sectors[i];
            }, py::return_value_policy::reference);

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

}