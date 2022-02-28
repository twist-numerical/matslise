#ifndef MATSLISE_JSMODULE
#define MATSLISE_JSMODULE

#include <emscripten/bind.h>
#include "../matslise.h"

using namespace emscripten;

inline val transformEigenvalues(const std::vector<std::pair<int, double>> &values) {
    val result{val::array()};

    for (std::pair<int, double> value: values)
        result.call<val>("push", val(value));

    return result;
}

inline Eigen::ArrayXd val2ArrayXd(const val &a) {
    int n = a[std::string("length")].as<int>();
    Eigen::ArrayXd r(n);
    for (int i = 0; i < n; ++i)
        r[i] = a[i].as<double>();
    return r;
}

template<typename T>
inline val vector2val(const std::vector<T> &a) {
    val rv = val::array();
    for (size_t i = 0; i < a.size(); ++i)
        rv.call<T>("push", a[i]);
    return rv;
}

inline val ArrayXd2val(const Eigen::ArrayXd &a) {
    val rv = val::array();
    for (Eigen::Index i = 0; i < a.size(); ++i)
        rv.call<double>("push", a[i]);
    return rv;
}

inline val ArrayXXd2val(const Eigen::ArrayXXd &a) {
    val r = val::array();
    for (Eigen::Index i = 0; i < a.rows(); ++i)
        r.call<val>("push", ArrayXd2val(a.row(i)));
    return r;
}

#endif // MATSLISE_JSMODULE