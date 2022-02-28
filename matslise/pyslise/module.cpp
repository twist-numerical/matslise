#include "module.h"

void pyslise(py::module &);

void pyscs(py::module &);
#ifndef WITH_MATSCS
inline void pyscs(py::module &) {};
#endif


PYBIND11_MODULE(pyslise, m) {
    pyslise(m);

    pyscs(m);
}
