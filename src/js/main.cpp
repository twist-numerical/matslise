#include <emscripten/bind.h>
#include "../matslise.h"
#include "../matslise2d.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

val transformEigenvalues(const vector<pair<int, double>> &values) {
    val result(val::array());

    for (pair<int, double> value: values)
        result.call<val>("push", val(value));

    return result;
}

ArrayXd val2ArrayXd(const val &a) {
    int n = a[string("length")].as<int>();
    ArrayXd r(n);
    for (int i = 0; i < n; ++i)
        r[i] = a[i].as<double>();
    return r;
}

template<typename T>
val vector2val(const vector<T> &a) {
    val rv = val::array();
    for (size_t i = 0; i < a.size(); ++i)
        rv.call<T>("push", a[i]);
    return rv;
}

val ArrayXd2val(const ArrayXd &a) {
    val rv = val::array();
    for (Eigen::Index i = 0; i < a.size(); ++i)
        rv.call<double>("push", a[i]);
    return rv;
}

val ArrayXXd2val(const ArrayXXd &a) {
    val r = val::array();
    for (Eigen::Index i = 0; i < a.rows(); ++i)
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

    class_<std::function<Y<>(double)>>("EigenfunctionCalculator")
            .function("eval",
                      optional_override([](std::function<Y<>(double)> &self, double x) -> Vector2d {
                          return self(x).y;
                      }));

    class_<std::function<double(double, double)>>("EigenfunctionCalculator2D")
            .function("eval",
                      optional_override([](std::function<double(double, double)> &self, double x, double y) -> double {
                          return self(x, y);
                      }));

    class_<AbstractMatslise<double>>("AbstractMatslise")
            .function("eigenfunction", optional_override(
                    [](const AbstractMatslise<double> &m, double E, const Vector2d &left, const Vector2d &right,
                       int index = -1) ->
                            val {
                        std::function<Y<>(double)> calculator =
                                m.eigenfunctionCalculator(E, Y<>(left, {0, 0}), Y<>(right, {0, 0}), index);
                        return val::global("Function")
                                .new_(string("calculator"), string(
                                        "var f = function(x) { return calculator.eval(x); };"
                                        "f.delete = function() { calculator.delete(); };"
                                        "return f;")
                                )(calculator);
                    }))
            .function("computeEigenfunction", optional_override(
                    [](const AbstractMatslise<double> &m, double E, const Vector2d &left, const Vector2d &right,
                       const val &x, int index = -1) ->
                            val {
                        Array<Y<double>, Dynamic, 1> array = m.eigenfunction(
                                E, Y<>(left, {0, 0}), Y<>(right, {0, 0}), val2ArrayXd(x), index);
                        return ArrayXd2val(array.unaryExpr([](const Y<double> &y) -> double {
                            return y.y(0);
                        }));
                    }))
            .function("eigenvaluesByIndex", optional_override(
                    [](const AbstractMatslise<double> &m, int Imin, int Imax, const Vector2d &left,
                       const Vector2d &right) -> val {
                        return transformEigenvalues(
                                m.eigenvaluesByIndex(Imin, Imax, Y<>(left, {0, 0}),
                                                     Y<>(right, {0, 0})));
                    }))
            .function("eigenvalueError", optional_override(
                    [](const AbstractMatslise<double> &m, double E, const Vector2d &left, const Vector2d &right,
                       int index = -1) -> double {
                        return m.eigenvalueError(E, Y<>(left, {0, 0}), Y<>(right, {0, 0}), index);
                    }));

    class_<Matslise<>, base<AbstractMatslise<double>>>("Matslise")
            .constructor(optional_override(
                    [](val f, double min, double max, const val &options) -> Matslise<> * {
                        std::shared_ptr<matslise::SectorBuilder<Matslise<>>> builder;
                        if (options["sectorCount"] != val::undefined())
                            builder = Matslise<>::UNIFORM(options["sectorCount"].as<int>());
                        else if (options["tolerance"] != val::undefined())
                            builder = Matslise<>::AUTO(options["tolerance"].as<double>());

                        return new Matslise<>([f](double x) -> double { return f(x).as<double>(); }, min, max, builder);
                    }))
            .function("propagate", optional_override(
                    [](const Matslise<> &m, double E, const Vector2d &y, double a, double b) ->
                            val {
                        Y<> y0;
                        double theta;
                        tie(y0, theta) = m.propagate(E, Y<>(y, {0, 0}), a, b);
                        val rv(val::object());
                        rv.set("y", val(y0.y));
                        rv.set("theta", val(theta));
                        return rv;
                    }))
            .function("sectorPoints", optional_override(
                    [](const Matslise<> &ms) -> val {
                        val r = val::array();
                        for (int i = 1; i < ms.sectorCount; ++i)
                            r.call<val>("push", ms.sectors[i]->min);
                        return r;
                    }));

    class_<MatsliseHalf<>, base<AbstractMatslise<double>>>("MatsliseHalfRange")
            .constructor(optional_override(
                    [](val f, double max, const val &options) -> MatsliseHalf<> * {
                        std::shared_ptr<matslise::SectorBuilder<Matslise<>>> builder;
                        if (options["sectorCount"] != val::undefined())
                            builder = Matslise<>::UNIFORM(options["sectorCount"].as<int>());
                        else if (options["tolerance"] != val::undefined())
                            builder = Matslise<>::AUTO(options["tolerance"].as<double>());

                        return new MatsliseHalf<>([f](double x) -> double { return f(x).as<double>(); }, max, builder);
                    }))
            .function("sectorPoints", optional_override(
                    [](const MatsliseHalf<> &m) -> val {
                        val r = val::array();
                        for (int i = 1; i < m.ms->sectorCount; ++i)
                            r.call<val>("push", m.ms->sectors[i]->min);
                        return r;
                    }));

    class_<Matslise2D<>>("Matslise2D")
            .constructor(optional_override(
                    [](val f, double xmin, double xmax, double ymin, double ymax, const val &options) -> Matslise2D<> * {
                        Options2<> o2;
                        if (options["sectorCount"] != val::undefined())
                            o2.sectorCount(options["sectorCount"].as<int>());
                        else if (options["tolerance"] != val::undefined())
                            o2.tolerance(options["tolerance"].as<double>());
                        Options1<> o1;
                        if (options["nested"] != val::undefined()) {
                            if (options["nested"]["sectorCount"] != val::undefined())
                                o1.sectorCount(options["nested"]["sectorCount"].as<int>());
                            else if (options["nested"]["tolerance"] != val::undefined())
                                o1.tolerance(options["nested"]["tolerance"].as<double>());
                        }
                        o2.nested(o1);
                        if (options["stepsPerSector"] != val::undefined())
                            o2.stepsPerSector(options["stepsPerSector"].as<int>());

                        return new Matslise2D<>([f](double x, double y) -> double { return f(x, y).as<double>(); },
                                                {{xmin, xmax}, ymin, ymax}, o2);
                    }))
            .function("eigenvaluesByIndex", optional_override([](const Matslise2D<> &se2d, int imin, int imax) -> val {
                return vector2val(se2d.eigenvaluesByIndex(imin, imax));
            }))
            .function("firstEigenvalue", &Matslise2D<double>::firstEigenvalue)
            .function("calculateError", optional_override([](Matslise2D<> &se2d, double E) -> pair<double, double> {
                return se2d.matchingError(E);
            }))
            .function("calculateErrors", optional_override([](Matslise2D<> &se2d, double E) -> val {
                vector<pair<double, double>> result = se2d.matchingErrors(E);
                val r = val::array();
                for (pair<double, double> &eigenvalue  :result)
                    r.call<val>("push", eigenvalue);
                return r;
            }))
            .function("eigenvalue", optional_override([](Matslise2D<> &se2d, double E) -> double {
                return se2d.eigenvalue(E);
            }))
            .function("sectorPoints", optional_override([](const Matslise2D<> &se2d) -> val {
                val r = val::array();
                for (int i = 1; i < se2d.sectorCount; ++i)
                    r.call<val>("push", se2d.sectors[i]->min);
                return r;
            }))
            .function("computeEigenfunction",
                      optional_override([](Matslise2D<> &se2d, double E, const val &x, const val &y) -> val {
                          vector<ArrayXXd> result = se2d.eigenfunction(E, val2ArrayXd(x), val2ArrayXd(y));
                          val r = val::array();
                          for (ArrayXXd &eigenfunction : result)
                              r.call<val>("push", ArrayXXd2val(eigenfunction));
                          return r;
                      }))
            .function("eigenfunction", optional_override(
                    [](const Matslise2D<> &se2d, double E) -> val {
                        std::vector<std::function<double(double, double)>> calculators = se2d.eigenfunctionCalculator(
                                E);
                        val r = val::array();
                        for (const auto &f : calculators)
                            r.call<val>("push", val::global("Function")
                                    .new_(string("calculator"), string(
                                            "var f = function(x, y) { return calculator.eval(x, y); };"
                                            "f.delete = function() { calculator.delete(); };"
                                            "return f;")
                                    )(f));
                        return r;
                    }));
}