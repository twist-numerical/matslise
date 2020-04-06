#include "module.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

void bind_matslise() {
    class_<std::function<Y<>(double)>>("EigenfunctionCalculator")
            .function("eval",
                      optional_override([](std::function<Y<>(double)> &self, double x) -> Vector2d {
                          return self(x).y;
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
                        matslise::SectorBuilder<Matslise<>> builder;
                        if (options["sectorCount"] != val::undefined())
                            builder = sector_builder::uniform<Matslise<>>(options["sectorCount"].as<int>());
                        else if (options["tolerance"] != val::undefined())
                            builder = sector_builder::automatic<Matslise<>>(options["tolerance"].as<double>());

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
                        matslise::SectorBuilder<Matslise<>> builder;
                        if (options["sectorCount"] != val::undefined())
                            builder = sector_builder::uniform<Matslise<>>(options["sectorCount"].as<int>());
                        else if (options["tolerance"] != val::undefined())
                            builder = sector_builder::automatic<Matslise<>>(options["tolerance"].as<double>());

                        return new MatsliseHalf<>([f](double x) -> double { return f(x).as<double>(); }, max, builder);
                    }))
            .function("sectorPoints", optional_override(
                    [](const MatsliseHalf<> &m) -> val {
                        val r = val::array();
                        for (int i = 1; i < m.ms->sectorCount; ++i)
                            r.call<val>("push", m.ms->sectors[i]->min);
                        return r;
                    }));
}