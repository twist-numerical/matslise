#include "module.h"

void pyslise(py::module &);

void pyslise_sturm_liouville(py::module &);

void pyscs(py::module &);


#ifndef WITH_MATSCS
inline void pyscs(py::module &) {};
#endif

#define xstr(s) as_str(s)
#define as_str(s) #s

PYBIND11_MODULE(pyslise, m) {
    m.attr("version") = py::str(xstr(MATSLISE_VERSION));

    pyslise(m);

    pyslise_sturm_liouville(m);

    pyscs(m);
}
