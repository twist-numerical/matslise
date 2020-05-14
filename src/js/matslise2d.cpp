#include "module.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

void bind_matslise2d() {
    class_<AbstractMatslise2D<double>>("AbstractMatslise2D")
            .function("firstEigenvalue", &AbstractMatslise2D<double>::firstEigenvalue)
            .function("eigenvalue", &AbstractMatslise2D<double>::eigenvalue)
            .function("eigenvalueError", &AbstractMatslise2D<double>::eigenvalueError)
            .function("eigenvalues", optional_override([](
                    const AbstractMatslise2D<double> &se2d, double emin, double emax) -> val {
                return vector2val(se2d.eigenvalues(emin, emax));
            }))
            .function("eigenvaluesByIndex",
                      optional_override([](const AbstractMatslise2D<double> &se2d, int imin, int imax) -> val {
                          return vector2val(se2d.eigenvaluesByIndex(imin, imax));
                      }))
            .function("computeEigenfunction",
                      optional_override(
                              [](const AbstractMatslise2D<double> &se2d, double E, const val &x, const val &y) -> val {
                                  vector<ArrayXXd> result = se2d.eigenfunction(E, val2ArrayXd(x), val2ArrayXd(y));
                                  val r = val::array();
                                  for (ArrayXXd &eigenfunction : result)
                                      r.call<val>("push", ArrayXXd2val(eigenfunction));
                                  return r;
                              }))
            .function("eigenfunction", optional_override(
                    [](const AbstractMatslise2D<double> &se2d, double E) -> val {
                        std::vector<std::function<double(double, double)>> calculators = se2d.eigenfunction(E);
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

    class_<Matslise2D<>, base<AbstractMatslise2D<double>>>("Matslise2D")
            .constructor(optional_override(
                    [](val f, double xmin, double xmax, double ymin, double ymax,
                       const val &options) -> Matslise2D<> * {
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

                            if (options["nested"]["symmetric"] != val::undefined())
                                o1.symmetric(options["nested"]["symmetric"].as<bool>());
                        }
                        o2.nested(o1);

                        return new Matslise2D<>([f](double x, double y) -> double { return f(x, y).as<double>(); },
                                                {{xmin, xmax}, ymin, ymax}, o2);
                    }))
            .function("sectorPoints", optional_override([](const Matslise2D<> &se2d) -> val {
                val r = val::array();
                for (int i = 1; i < se2d.sectorCount; ++i)
                    r.call<val>("push", se2d.sectors[i]->min);
                return r;
            }))
            .function("matchingError", optional_override([](Matslise2D<> &se2d, double E) -> pair<double, double> {
                return se2d.matchingError(E);
            }))
            .function("matchingErrors", optional_override([](Matslise2D<> &se2d, double E) -> val {
                vector<pair<double, double>> result = se2d.matchingErrors(E);
                val r = val::array();
                for (pair<double, double> &eigenvalue  :result)
                    r.call<val>("push", eigenvalue);
                return r;
            }));

    class_<Matslise2DHalf<>>("Matslise2DHalf")
            .constructor(optional_override(
                    [](val f, double xmin, double xmax, double ymin, double ymax,
                       const val &options) -> Matslise2DHalf<> * {
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

                            if (options["nested"]["symmetric"] != val::undefined())
                                o1.symmetric(options["nested"]["symmetric"].as<bool>());
                        }
                        o2.nested(o1);

                        return new Matslise2DHalf<>([f](double x, double y) -> double { return f(x, y).as<double>(); },
                                                    {{xmin, xmax}, ymin, ymax}, o2);
                    }));
}