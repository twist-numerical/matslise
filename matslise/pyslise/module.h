#ifndef MATSLISE_PYSLISE
#define MATSLISE_PYSLISE

#include "../matslise.h"
#include "../matscs.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include <tuple>

namespace py = pybind11;
using namespace Eigen;
using namespace std;
using namespace matslise;
using namespace sector_builder;

template<int cols = 1>
inline Y<double, 1, cols> make_y(const Matrix<double, 2, cols, 0, 2, cols> &value) {
    return Y<double, 1, cols>(value, Matrix<double, 2, cols>::Zero());
}

inline Y<> packY(const tuple<double, double> &t) {
    return make_y<1>({get<0>(t), get<1>(t)});
}

inline Y<double, Dynamic> packY(const tuple<MatrixXd, MatrixXd> &t) {
    MatrixXd a, b;
    tie(a, b) = t;
    Y<double, Dynamic> result(a.rows(), a.cols());
    result.block() = a;
    result.block(dX) = b;
    return result;
}

inline Y<> packY(const tuple<double, double> &t, const tuple<double, double> &dt) {
    return Y<>({get<0>(t), get<1>(t)},
               {get<0>(dt), get<1>(dt)});
}

inline pair<pair<MatrixXd, MatrixXd>, pair<MatrixXd, MatrixXd>> unpackY(const Y<double, Dynamic> &y) {
    return make_pair(make_pair(y.block(), y.block(dX)), make_pair(y.block(dE), y.block(dXdE)));
}

inline pair<pair<double, double>, pair<double, double>> unpackY(const Y<> &y) {
    return make_pair(make_pair(y.y()[0], y.y()[1]), make_pair(y.ydE()[0], y.ydE()[1]));
}

template<typename Scalar, typename Nurse>
class EigenfunctionWrapper : public AbstractMatslise<Scalar>::Eigenfunction {
public:
    std::shared_ptr<Nurse> nurse;
    std::unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction> function;

    EigenfunctionWrapper(
            std::shared_ptr<Nurse> nurse_,
            std::unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction> function_) :
            nurse(std::move(nurse_)), function(std::move(function_)) {
    }

    virtual Eigen::Array<Scalar, 2, 1> operator()(const Scalar &x) const override {
        return (*function)(x);
    };

    virtual Eigen::Array<Scalar, Eigen::Dynamic, 2>
    operator()(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const override {
        return (*function)(x);
    }
};

#endif
