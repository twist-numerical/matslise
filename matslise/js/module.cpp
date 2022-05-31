#include "module.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

void bind_matslise();

void bind_matslise2d();

EMSCRIPTEN_BINDINGS(module) {
    value_array<Vector2d>("Vector2d")
            .element(emscripten::index<0>())
            .element(emscripten::index<1>());

    value_object<Eigenvalue>("Eigenvalue")
            .field("index", &Eigenvalue::index)
            .field("eigenvalue", &Eigenvalue::eigenvalue);

    value_object<Eigenpair>("Eigenpair")
            .field("index", &Eigenpair::index)
            .field("eigenvalue", &Eigenpair::eigenvalue)
            .field("eigenfunction", &Eigenpair::eigenfunction);

    value_object<pair<double, double>>("PairDoubleDouble")
            .field("first", &pair<double, double>::first)
            .field("second", &pair<double, double>::second);

    bind_matslise();
}
