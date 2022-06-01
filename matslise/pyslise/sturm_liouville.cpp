#include <utility>

#include "../liouville.h"
#include "module.h"

void pyslise_sturm_liouville(py::module &m) {
    py::class_<LiouvilleTransformation<double>>(m, "LiouvilleTransformation")
            .def(py::init([](const function<double(double)> &p, const function<double(double)> &q,
                             const function<double(double)> &w, double min, double max) {
                return LiouvilleTransformation<double>({min, max}, p, q, w);
            }), R""""(\
)"""", py::arg("p"), py::arg("q"), py::arg("w"), py::arg("min"), py::arg("max"))
            .def("V", &LiouvilleTransformation<double>::V, py::arg("x"))
            .def("r2x", &LiouvilleTransformation<double>::r2x, py::arg("r"))
            .def("x2r", &LiouvilleTransformation<double>::x2r, py::arg("x"))
            .def_property_readonly("rDomain", [](const LiouvilleTransformation<double> &lt) {
                return std::pair{lt.rDomain().min(), lt.rDomain().max()};
            })
            .def_property_readonly("xDomain", [](const LiouvilleTransformation<double> &lt) {
                auto domain = lt.xDomain();
                return std::pair{domain.min(), domain.max()};
            })
            .def_property_readonly("pieces", [](const LiouvilleTransformation<double> &lt) {
                std::vector<std::pair<double, double>> r;
                r.reserve(lt.pieces.size());
                std::transform(lt.pieces.begin(), lt.pieces.end(), std::back_inserter(r), [](auto p) {
                    return std::pair{p.r.min(), p.r.max()};
                });
                return r;
            });

    py::class_<SturmLiouville<double>, shared_ptr<SturmLiouville<double>>>(m, "SturmLiouville", R""""(\
>>> from math import pi, sqrt, sin
>>> import numpy as np
>>> slp = SturmLiouville(lambda x: 2+sin(2*pi*x), lambda x: -10, lambda x: 1 + sqrt(x), 0, 1, 1e-8)
>>> i, E = slp.eigenvaluesByIndex(2, 3, (0,1), (2, -10))[0]
>>> abs(E - 66.9259933904942) < 1e-6
True
>>> i, E, f = slp.eigenpairsByIndex(1, 2, (0, 1), (2, -10))[0]
>>> abs(E - 24.0804524555819) < 1e-6
True
>>> abs(f(0.27)[0] - 0.82637509455759) < 1e-6
True
)"""")
            .def(py::init([](const function<double(double)> &p, const function<double(double)> &q,
                             const function<double(double)> &w, double min, double max, double tolerance) {
                     return std::make_shared<SturmLiouville<double>>(p, q, w, Rectangle<double, 1>{min, max}, tolerance);
                 }), R""""()"""", py::arg("p"), py::arg("q"), py::arg("w"), py::arg("min"), py::arg("max"),
                 py::arg("tolerance") = 1e-8)
            .def("eigenvaluesByIndex",
                 [](const shared_ptr<SturmLiouville<double>> &slp, int Imin, int Imax, const Vector2d &left,
                    const optional<Vector2d> &_right) -> vector<pair<int, double>> {
                     const Vector2d &right = _right ? *_right : left;
                     return slp->eigenvaluesByIndex(Imin, Imax, make_y(left), make_y(right));
                 }, R""""()"""",
                 py::arg("Imin"), py::arg("Imax"), py::arg("left"), py::arg("right") = optional<Vector2d>())
            .def("eigenpairsByIndex",
                 [](const shared_ptr<SturmLiouville<double>> &slp, int iMin, int iMax, const Vector2d &left,
                    const optional<Vector2d> &_right) {
                     const Vector2d &right = _right ? *_right : left;
                     vector<tuple<int, double, unique_ptr<AbstractMatslise<double>::Eigenfunction>>> result;
                     auto pairs = slp->eigenpairsByIndex(iMin, iMax, make_y(left), make_y(right));
                     result.reserve(pairs.size());
                     for (auto &eigenpair: pairs) {
                         result.emplace_back(get<0>(eigenpair), get<1>(eigenpair),
                                             std::make_unique<EigenfunctionWrapper<double, SturmLiouville<double>>>
                                                     (slp, std::move(get<2>(eigenpair))));
                     }
                     return result;
                 }, py::arg("Imin"), py::arg("Imax"), py::arg("left"), py::arg("right") = optional<Vector2d>());
}
