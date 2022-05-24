#include "module.h"

void pyslise(py::module &);

void pyslise_sturm_liouville(py::module &);

void pyscs(py::module &);


#ifndef WITH_MATSCS
inline void pyscs(py::module &) {};
#endif


PYBIND11_MODULE(pyslise, m) {
    pyslise(m);

    pyslise_sturm_liouville(m);

    pyscs(m);
}
