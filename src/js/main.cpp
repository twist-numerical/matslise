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
            .field("index", &pair<int, double>::first)
            .field("value", &pair<int, double>::second);

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
                        tie(y0, theta) = m.propagate(E, Y<double>(y), a, b);
                        val rv(val::object());
                        rv.set("y", val(y0.y.toEigen()));
                        rv.set("theta", val(theta));
                        return rv;
                    }))
            .function("eigenfunction", optional_override(
                    [](Matslise &m, double E, const Vector2d &left, const Vector2d &right) ->
                            val {
                        matslise_util::EigenfunctionCalculator *calculator =
                                m.eigenfunctionCalculator(E, Y<double>(left), Y<double>(right));
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
                                m.computeEigenvaluesByIndex(Imin, Imax, Y<double>(left), Y<double>(right)));
                    }), allow_raw_pointers());

    class_<se2d_util::Sector>("Sector2D")
            .function("computeEigenfunction",
                      optional_override([](const se2d_util::Sector &sector, int index, val x) -> val {
                          return ArrayXd2val(sector.computeEigenfunction(index, val2ArrayXd(x)));
                      }))
            .function("computeEigenfunctions",
                      optional_override([](const se2d_util::Sector &sector, int index, val x) -> val {
                          ArrayXd ax = val2ArrayXd(x);
                          val r = val::array();
                          for (int i = 0; i < sector.se2d->N; ++i)
                              r.call<val>("push", ArrayXd2val(sector.computeEigenfunction(index, ax)));
                          return r;
                      }));

    class_<SE2D>("SE2D")
            .constructor(optional_override(
                    [](val f, double xmin, double xmax, double ymin, double ymax, int xSectorCount,
                       int ySectorCount) -> SE2D * {
                        return new SE2D([f](double x, double y) -> double { return f(x, y).as<double>(); }, xmin, xmax,
                                        ymin, ymax, xSectorCount, ySectorCount);
                    }))
            .function("calculateError", &SE2D::calculateError)
            .function("computeEigenfunction", optional_override([](SE2D &se2d, double E, val x, val y) -> val {
                vector<ArrayXXd> result = se2d.computeEigenfunction(E, val2ArrayXd(x), val2ArrayXd(y));
                val r = val::array();
                for (ArrayXXd &eigenfunction  : result)
                    r.call<val>("push", ArrayXXd2val(eigenfunction));
                return r;
            }))
            .function("computeEigenfunctionSteps", optional_override([](SE2D &se2d, double E) -> val {
                Y<MatrixXd> *result = se2d.computeEigenfunctionSteps(E);
                val r = val::array();
                for (int i = 0; i <= se2d.sectorCount; ++i) {
                    val y = val::array();
                    y.call<val>("push", ArrayXXd2val(result[i].y.x.array()));
                    y.call<val>("push", ArrayXXd2val(result[i].y.y.array()));
                    r.call<val>("push", y);
                }
                return r;
            }))
            .function("getSector", &SE2D::getSector);
}