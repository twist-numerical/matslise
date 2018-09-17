#include <emscripten/bind.h>
#include "../matslise.h"
#include "../se2d.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

val transformEigenvalues(vector<pair<unsigned int, double>> *values) {
    val result(val::array());

    for (pair<unsigned int, double> value: *values)
        result.call<val>("push", val(value));

    delete values;
    return result;
}

EMSCRIPTEN_BINDINGS(Matslise) {
    value_array<Vector2d>("Vector2d")
            .element(emscripten::index<0>())
            .element(emscripten::index<1>());

    value_object<pair<unsigned int, double>>("PairIntDouble")
            .field("index", &pair<unsigned int, double>::first)
            .field("value", &pair<unsigned int, double>::second);

    class_<matslise_util::EigenfunctionCalculator>("EigenfunctionCalculator")
            .function("eval",
                      optional_override([](matslise_util::EigenfunctionCalculator &self, double x) -> Vector2d { return self(x).y; }));

    class_<Matslise>("Matslise")
            .constructor(optional_override(
                    [](val f, double min, double max, int steps) -> Matslise * {
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
                    }))
            .function("eigenfunction", optional_override(
                    [](Matslise &m, double E, const Vector2d &left, const Vector2d &right) ->
                            val {
                        matslise_util::EigenfunctionCalculator *calculator = m.eigenfunctionCalculator(E, Y<double>(left),
                                                                                        Y<double>(right));
                        return val::global("Function")
                                .new_(string("calculator"), string(
                                        "var f = function(x) { return calculator.eval(x); };"
                                        "f.delete = function() { calculator.delete(); };"
                                        "return f;")
                                )(calculator);
                    }))
            .function("eigenvaluesByIndex", optional_override(
                    [](Matslise &m, unsigned int Imin, unsigned int Imax, const Vector2d &left,
                       const Vector2d &right) -> val {
                        return transformEigenvalues(
                                m.computeEigenvaluesByIndex(Imin, Imax, Y<double>(left), Y<double>(right)));
                    }), allow_raw_pointers());
}