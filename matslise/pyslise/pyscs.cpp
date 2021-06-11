#include "module.h"

void pyscs(py::module &m) {
    py::class_<Matscs<>::Sector>(m, "PyscsSector")
            .def_readonly("min", &Matscs<>::Sector::min)
            .def_readonly("max", &Matscs<>::Sector::max)
            .def_property_readonly("forward", [](const Matscs<>::Sector &s) -> bool {
                return s.direction == Direction::forward;
            });
    py::class_<Matscs<>>(m, "Pyscs")
            .def(py::init(
                    [](const function<MatrixXd(double)> &V, int N, double xmin, double xmax, int steps,
                       double tolerance) {
                        if (steps != -1 && tolerance != -1)
                            throw invalid_argument("Not both 'steps' and 'tolerance' can be set.");
                        if (steps == -1 && tolerance == -1)
                            throw invalid_argument("One of 'steps' and 'tolerance' must be set.");
                        return new Matscs<>(
                                V, N, xmin, xmax,
                                steps != -1 ? sector_builder::uniform<Matscs<double>>(steps)
                                            : sector_builder::automatic<Matscs<double>>(tolerance));
                    }), "Pyscs", py::arg("V"), py::arg("dimensions"), py::arg("xmin"), py::arg("xmax"),
                 py::arg("steps") = -1,
                 py::arg("tolerance") = -1)
            .def_readonly("__sectorCount", &Matscs<>::sectorCount)
            .def_readonly("__matchIndex", &Matscs<>::matchIndex)
            .def_readonly("min", &Matscs<>::xmin)
            .def_readonly("max", &Matscs<>::xmax)
            .def("propagate",
                 [](Matscs<> &m, double E, const tuple<MatrixXd, MatrixXd> &y, double a,
                    double b) -> pair<pair<MatrixXd, MatrixXd>, double> {
                     double theta = 0;
                     Y<double, Dynamic> r;
                     tie(r, theta) = m.propagate(E, packY(y), a, b, true);
                     return make_pair(make_pair(r.block(), r.block(dX)), theta);
                 })
            .def("propagatePsi", &Matscs<>::propagatePsi)
            .def("__sector", [](Matscs<> &p, int i) -> Matscs<>::Sector * {
                return p.sectors[i].get();
            }, py::return_value_policy::reference);
}