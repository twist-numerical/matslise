#ifndef MATSLISE_JSMODULE
#define MATSLISE_JSMODULE

#include <emscripten/bind.h>
#include "../matslise.h"

using namespace emscripten;

struct Eigenvalue {
    int index = -1;
    double eigenvalue = 0;

    Eigenvalue() = default;

    Eigenvalue(int index_, double eigenvalue_) : index(index_), eigenvalue(eigenvalue_) {}
};

struct Eigenpair : public Eigenvalue {
    val eigenfunction = val::undefined();

    Eigenpair() = default;

    Eigenpair(int index_, double eigenvalue_, val eigenfunction_)
            : Eigenvalue(index_, eigenvalue_), eigenfunction(eigenfunction_) {}
};

template<typename It, typename F>
inline val toValArray(It begin, It end, F f = [](auto x) { return x; }) {
    val result{val::array()};

    std::for_each(begin, end, [&f, &result](auto &t) {
        result.call<val>("push", val(f(t)));
    });

    return result;
}

inline Eigen::ArrayXd val2ArrayXd(const val &a) {
    int n = a[std::string("length")].as<int>();
    Eigen::ArrayXd r(n);
    for (int i = 0; i < n; ++i)
        r[i] = a[i].as<double>();
    return r;
}

inline val ArrayXd2val(const Eigen::ArrayXd &a) {
    val rv = val::array();
    for (Eigen::Index i = 0; i < a.size(); ++i)
        rv.call<double>("push", a[i]);
    return rv;
}

inline val ArrayX2d2val(const Eigen::ArrayX2d &a) {
    val rv = val::array();
    for (Eigen::Index i = 0; i < a.rows(); ++i)
        rv.call<double>("push", (Eigen::Vector2d) a.row(i));
    return rv;
}

#endif // MATSLISE_JSMODULE
