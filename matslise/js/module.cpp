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

    value_object<pair<int, double>>("PairIntDouble")
            .field("first", &pair<int, double>::first)
            .field("second", &pair<int, double>::second);

    value_object<pair<double, double>>("PairDoubleDouble")
            .field("first", &pair<double, double>::first)
            .field("second", &pair<double, double>::second);

    bind_matslise();
}
