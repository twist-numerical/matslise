#include <emscripten/bind.h>
#include "../matslise.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;


EMSCRIPTEN_BINDINGS(Matslise) {
    value_array<Vector2d>("Vector2d")
            .element(emscripten::index<0>())
            .element(emscripten::index<1>());

    class_<Matslise>("Matslise")
            .constructor(optional_override(
                    [](val f, double min, double max, int steps) {
                        return new Matslise([f](double x) -> double { return f(x).as<double>(); }, min, max, steps);
                    }))
            .function("propagate", optional_override(
                    [](Matslise &m, double E, const Vector2d &y, double a, double b) ->
                            val {
                        Y<double> y0;
                        double theta;
                        tie(y0, theta) = m.propagate(E, Y<double>(y), a, b);
                        val rv(val::object());
                        rv.set("y", val((Vector2d) (y0.y)));
                        rv.set("theta", val(theta));
                        return rv;
                    }));
}