#include <pybind11/pybind11.h>
#include <matslise.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

namespace py = pybind11;

matslise::Y packY(std::tuple<double, double> t) {
    return matslise::Y(std::get<0>(t), std::get<1>(t));
}
std::tuple<double, double> unpackY(matslise::Y y) {
    std::tuple<double, double> t(y.y, y.dy);
    return t;
};

PYBIND11_MODULE(pyslise, m) {
    py::class_<Matslise>(m, "Pyslise")
            .def(py::init<std::function<double(double)>, double, double, int>())
            .def("propagate",
                [](Matslise &m, double E, std::tuple<double, double> y, double a, double b) ->
                    std::tuple<double, double> {
                        return unpackY(m.propagate(E, packY(y), a, b));
                    })
            .def("computeEigenfunction", [](Matslise &m, double E, std::vector<double> xs)
                 -> std::tuple<std::vector<double>, std::vector<std::tuple<double, double>>*> {
                    auto ysY = m.computeEigenfunction(E, xs);
                    auto ys = new std::vector<std::tuple<double, double>>();
                    for(matslise::Y y : *ysY)
                        ys->push_back(unpackY(y));
                    delete ysY;
                    std::tuple<std::vector<double>, std::vector<std::tuple<double, double>>*> xsys(xs, ys);
                    return xsys;
                });

}