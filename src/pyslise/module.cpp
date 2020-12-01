#include "module.h"

void pyslise(py::module &);

void pyscs(py::module &);

void pyslise2d(py::module &);

void pyslise3d(py::module &);

PYBIND11_MODULE(pyslise, m) {
    pyslise(m);

    pyscs(m);

    pyslise2d(m);

    pyslise3d(m);
}
