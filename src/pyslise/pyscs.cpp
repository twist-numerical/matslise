#include "module.h"

void pyScs(py::module &m) {
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
                    double b) -> pair<pair<MatrixXd, MatrixXd>, double> {
                     double theta = 0;
                     Y<double, Dynamic> r;
                     tie(r, theta) = m.propagate(E, packY(y), a, b, true);
                     return make_pair(make_pair(r.getY(0), r.getY(1)), theta);
                 })
            .def("propagatePsi", &Matscs<>::propagatePsi)
            .def("__sector", [](Matscs<> &p, int i) -> Matscs<>::Sector * {
                return p.sectors[i];
            }, py::return_value_policy::reference);
}