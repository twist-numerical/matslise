#include "module.h"

void pyslise(py::module &);

void pyscs(py::module &);
#ifndef WITH_MATSCS
inline void pyscs(py::module &) {};
#endif

void pyslise2d(py::module &);
#ifndef WITH_MATSLISE_2D
inline void pyslise2d(py::module &) {};
#endif

void pyslise3d(py::module &);
#ifndef WITH_MATSLIS_3D
inline void pyslise3d(py::module &) {};
#endif

PYBIND11_MODULE(pyslise, m) {
    pyslise(m);

    pyscs(m);

    pyslise2d(m);

    pyslise3d(m);
}
