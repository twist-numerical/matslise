#include <emscripten/bind.h>
#include "../matslise.h"
#include "../se2d.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

val transformEigenvalues(vector<pair<int, double>> *values) {
    val result(val::array());

    for (pair<int, double> value: *values)
        result.call<val>("push", val(value));

    delete values;
    return result;
}

ArrayXd val2ArrayXd(const val &a) {
    int n = a[string("length")].as<int>();
    ArrayXd r(n);
    for (int i = 0; i < n; ++i)
        r[i] = a[i].as<double>();
    return r;
}

val ArrayXd2val(const ArrayXd &a) {
    val rv = val::array();
    for (int i = 0; i < a.size(); ++i)
        rv.call<double>("push", a[i]);
    return rv;
}

val ArrayXXd2val(const ArrayXXd &a) {
    val r = val::array();
    for (int i = 0; i < a.rows(); ++i)
        r.call<val>("push", ArrayXd2val(a.row(i)));
    return r;
}

EMSCRIPTEN_BINDINGS(Matslise) {
    value_array<Vector2d>("Vector2d")
            .element(emscripten::index<0>())
            .element(emscripten::index<1>());

    value_object<pair<int, double>>("PairIntDouble")
            .field("first", &pair<int, double>::first)
            .field("second", &pair<int, double>::second);

    value_object<pair<double, double>>("PairDoubleDouble")
            .field("first", &pair<double, double>::first)
            .field("second", &pair<double, double>::second);

    class_<matslise_util::EigenfunctionCalculator>("EigenfunctionCalculator")
            .function("eval",
                      optional_override([](matslise_util::EigenfunctionCalculator &self, double x) -> Vector2d {
                          return self(x).y.toEigen();
                      }));

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
                        tie(y0, theta) = m.propagate(E, Y<double>(y, {0, 0}), a, b);
                        val rv(val::object());
                        rv.set("y", val(y0.y.toEigen()));
                        rv.set("theta", val(theta));
                        return rv;
                    }))
            .function("eigenfunction", optional_override(
                    [](Matslise &m, double E, const Vector2d &left, const Vector2d &right) ->
                            val {
                        matslise_util::EigenfunctionCalculator *calculator =
                                m.eigenfunctionCalculator(E, Y<double>(left, {0, 0}), Y<double>(right, {0, 0}));
                        return val::global("Function")
                                .new_(string("calculator"), string(
                                        "var f = function(x) { return calculator.eval(x); };"
                                        "f.delete = function() { calculator.delete(); };"
                                        "return f;")
                                )(calculator);
                    }))
            .function("eigenvaluesByIndex", optional_override(
                    [](Matslise &m, int Imin, int Imax, const Vector2d &left,
                       const Vector2d &right) -> val {
                        return transformEigenvalues(
                                m.computeEigenvaluesByIndex(Imin, Imax, Y<double>(left, {0, 0}),
                                                            Y<double>(right, {0, 0})));
                    }), allow_raw_pointers());

    class_<SEnD<2>>("SE2D")
            .constructor(optional_override(
                    [](val f, double xmin, double xmax, double ymin, double ymax,
                       int xSectorCount, int ySectorCount) -> SEnD<2> * {
                        return new SEnD<2>([f](double x, double y) -> double { return f(x, y).as<double>(); },
                                           {{xmin, xmax}, ymin, ymax},
                                           Options<2>().sectorCount(ySectorCount).nested(
                                                   Options<1>().sectorCount(xSectorCount)));
                    }))
            .function("calculateError", optional_override([](SEnD<2> &se2d, double E) -> std::pair<double, double> {
                return se2d.calculateError(E, SEnD_util::ABS_SORTER);
            }))
            .function("findEigenvalue", optional_override([](SEnD<2> &se2d, double E) -> double {
                return se2d.findEigenvalue(E);
            }))
            .function("computeEigenfunction", optional_override([](SEnD<2> &se2d, double E, val x, val y) -> val {
                vector<ArrayXXd> *result = se2d.computeEigenfunction(E, val2ArrayXd(x), val2ArrayXd(y));
                val r = val::array();
                for (ArrayXXd &eigenfunction  : *result)
                    r.call<val>("push", ArrayXXd2val(eigenfunction));
                delete result;
                return r;
            }));
}