#include "module.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

val eigenvaluesToVal(const vector<tuple<Index, double, Index>> &list) {
    val array = val::array();
    for (auto &item : list) {
        val objItem = val::object();
        objItem.set("index", get<0>(item));
        objItem.set("value", get<1>(item));
        objItem.set("multiplicity", get<2>(item));
        array.call<void>("push", objItem);
    }
    return array;
}

void bind_matslise2d() {
    class_<AbstractMatslise2D<double>>("AbstractMatslise2D")
            .function("eigenvalue", &AbstractMatslise2D<double>::eigenvalue)
            .function("eigenvalueError", &AbstractMatslise2D<double>::eigenvalueError)
            .function("eigenvalues", optional_override([](
                    const AbstractMatslise2D<double> &se2d, double emin, double emax) -> val {
                return eigenvaluesToVal(se2d.eigenvalues(emin, emax));
            }))
            .function("eigenvaluesByIndex",
                      optional_override([](const AbstractMatslise2D<double> &se2d, int imin, int imax) -> val {
                          return eigenvaluesToVal(se2d.eigenvaluesByIndex(imin, imax));
                      }))
            .function("computeEigenfunction",
                      optional_override(
                              [](const AbstractMatslise2D<double> &se2d, double E, const val &_x,
                                 const val &_y) -> val {
                                  ArrayXd x = val2ArrayXd(_x);
                                  ArrayXd y = val2ArrayXd(_y);
                                  val result = val::array();
                                  for (const auto &f : se2d.eigenfunction(E))
                                      result.call<val>("push", ArrayXXd2val(f(x, y)));
                                  return result;
                              }))
            .function("eigenfunction", optional_override(
                    [](const AbstractMatslise2D<double> &se2d, double E) -> val {
                        std::vector<Eigenfunction2D<>> calculators = se2d.eigenfunction(E);
                        val r = val::array();
                        for (const auto &f : calculators)
                            r.call<val>("push", val::global("Function")
                                    .new_(string("calculator"), string(
                                            "var f = function(x, y) { return calculator.eval(x, y); };"
                                            "f.delete = function() { calculator.delete(); };"
                                            "return f;")
                                    )((std::function<double(double, double)>) f));
                        return r;
                    }));

    class_<Matslise2D<>, base<AbstractMatslise2D<double>>>("Matslise2D")
            .constructor(optional_override(
                    [](val f, double xmin, double xmax, double ymin, double ymax,
                       const val &options) -> Matslise2D<> * {
                        Matslise2D<>::Config config;
                        if (options["tolerance"] != val::undefined())
                            config.tolerance = options["tolerance"].as<double>();

                        if (options["xSymmetric"] != val::undefined())
                            config.xSymmetric = options["xSymmetric"].as<bool>();

                        if (options["xSectorCount"] != val::undefined())
                            config.xSectorBuilder = sector_builder::uniform<Matslise<>>(
                                    options["xSectorCount"].as<int>());
                        else if (options["xTolerance"] != val::undefined())
                            config.xSectorBuilder = sector_builder::automatic<Matslise<>>(
                                    options["xTolerance"].as<double>());

                        if (options["ySectorCount"] != val::undefined())
                            config.ySectorBuilder = sector_builder::uniform<Matslise2D<>>(
                                    options["ySectorCount"].as<int>());
                        else if (options["yTolerance"] != val::undefined())
                            config.ySectorBuilder = sector_builder::automatic<Matslise2D<>>(
                                    options["yTolerance"].as<double>());

                        if (options["xSectorCount"] != val::undefined())
                            config.xSectorBuilder = sector_builder::uniform<Matslise<>>(
                                    options["xSectorCount"].as<int>());

                        return new Matslise2D<>([f](double x, double y) -> double { return f(x, y).as<double>(); },
                                                {{xmin, xmax}, ymin, ymax}, config);
                    }))
            .function("sectorPoints", optional_override([](const Matslise2D<> &se2d) -> val {
                val r = val::array();
                for (unsigned long i = 1; i < se2d.sectors.size(); ++i)
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
                        Matslise2D<>::Config config;
                        if (options["tolerance"] != val::undefined())
                            config.tolerance = options["tolerance"].as<double>();

                        if (options["xSymmetric"] != val::undefined())
                            config.xSymmetric = options["xSymmetric"].as<bool>();

                        if (options["xSectorCount"] != val::undefined())
                            config.xSectorBuilder = sector_builder::uniform<Matslise<>>(
                                    options["xSectorCount"].as<int>());
                        else if (options["xTolerance"] != val::undefined())
                            config.xSectorBuilder = sector_builder::automatic<Matslise<>>(
                                    options["xTolerance"].as<double>());

                        if (options["ySectorCount"] != val::undefined())
                            config.ySectorBuilder = sector_builder::uniform<Matslise2D<>>(
                                    options["ySectorCount"].as<int>());
                        else if (options["yTolerance"] != val::undefined())
                            config.ySectorBuilder = sector_builder::automatic<Matslise2D<>>(
                                    options["yTolerance"].as<double>());

                        if (options["xSectorCount"] != val::undefined())
                            config.xSectorBuilder = sector_builder::uniform<Matslise<>>(
                                    options["xSectorCount"].as<int>());

                        return new Matslise2DHalf<>([f](double x, double y) -> double { return f(x, y).as<double>(); },
                                                    {{xmin, xmax}, ymin, ymax}, config);
                    }));
}