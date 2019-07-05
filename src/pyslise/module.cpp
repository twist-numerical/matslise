#include "module.h"

void pySlise(py::module &);

void pyScs(py::module &);

void pySE2d(py::module &);

PYBIND11_MODULE(pyslise, m) {
    pySlise(m);

    pyScs(m);

    pySE2d(m);
}