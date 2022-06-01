#include "module.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

void bind_matslise() {

    class_<AbstractMatslise<double>::Eigenfunction>("Eigenfunction")
            .function("eval",
                      optional_override([](const AbstractMatslise<double>::Eigenfunction &self, val x) -> val {
                          if (val::global("Array").call<bool>("isArray", x)) {
                              return ArrayX2d2val(self(val2ArrayXd(x)));
                          } else {
                              return val((Vector2d) self(x.as<double>()));
                          }
                      }));

    class_<AbstractMatslise<double>>("AbstractMatslise")
            .function("eigenfunction", optional_override(
                    [](const AbstractMatslise<double> &m, double E, const Vector2d &left, const Vector2d &right,
                       int index = -1) -> val {
                        return wrapEigenfunction(m.eigenfunction(E, Y<>(left, {0, 0}), Y<>(right, {0, 0}), index));
                    }))
            .function("computeEigenfunction", optional_override(
                    [](const AbstractMatslise<double> &m, double E, const Vector2d &left, const Vector2d &right,
                       const val &x, int index = -1) -> val {
                        Array<double, Dynamic, 2> array = (*m.eigenfunction(
                                E, Y<>(left, {0, 0}), Y<>(right, {0, 0}), index))(val2ArrayXd(x));
                        return ArrayXd2val(array.col(0));
                    }))
            .function("eigenvaluesByIndex", optional_override(
                    [](const AbstractMatslise<double> &m, int Imin, int Imax, const Vector2d &left,
                       const Vector2d &right) -> val {
                        auto es = m.eigenvaluesByIndex(Imin, Imax, Y<>(left, {0, 0}), Y<>(right, {0, 0}));
                        return toValArray(es.begin(), es.end(), [](auto t) { return Eigenvalue(t.first, t.second); });
                    }))
            .function("eigenpairsByIndex", optional_override(
                    [](const AbstractMatslise<double> &m, int Imin, int Imax, const Vector2d &left,
                       const Vector2d &right) -> val {
                        auto es = m.eigenpairsByIndex(Imin, Imax, Y<>(left, {0, 0}), Y<>(right, {0, 0}));
                        return toValArray(es.begin(), es.end(), [](auto &t) {
                            return Eigenpair(get<0>(t), get<1>(t), wrapEigenfunction(std::move(get<2>(t))));
                        });
                    }))
            .function("eigenvalueError", optional_override(
                    [](const AbstractMatslise<double> &m, double E, const Vector2d &left, const Vector2d &right,
                       int index = -1) -> double {
                        return m.eigenvalueError(E, Y<>(left, {0, 0}), Y<>(right, {0, 0}), index);
                    }));

    class_<Matslise<>, base<AbstractMatslise<double>>>("Matslise")
            .constructor(optional_override(
                    [](val f, double min, double max, double tolerance) -> Matslise<> * {
                        return new Matslise<>(val2function<double(double)>(f), min, max, tolerance);
                    }))
            .function("propagate", optional_override(
                    [](const Matslise<> &m, double E, const Vector2d &y, double a, double b) ->
                            val {
                        Y<> y0;
                        double theta;
                        tie(y0, theta) = m.propagate(E, Y<>(y, {0, 0}), a, b);
                        val rv(val::object());
                        rv.set("y", val((Vector2d) y0.y()));
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

    class_<MatsliseHalf<>, base<AbstractMatslise<double>>>("MatsliseHalf")
            .constructor(optional_override(
                    [](val f, double max, double tolerance) -> MatsliseHalf<> * {
                        return new MatsliseHalf<>(val2function<double(double)>(f), max, tolerance);
                    }))
            .function("sectorPoints", optional_override(
                    [](const MatsliseHalf<> &m) -> val {
                        val r = val::array();
                        for (int i = 1; i < m.ms->sectorCount; ++i)
                            r.call<val>("push", m.ms->sectors[i]->min);
                        return r;
                    }));
}