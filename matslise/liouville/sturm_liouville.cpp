#include <vector>
#include "../liouville.h"

using namespace matslise;
using std::get;
using std::vector;
using std::unique_ptr;
using std::tuple;

template<typename Scalar>
class SLEigenfunction : public SturmLiouville<Scalar>::Eigenfunction {
public:
    const LiouvilleTransformation<Scalar> *transformation;
    const unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction> mEigenfunction;

public:
    SLEigenfunction(const LiouvilleTransformation<Scalar> *transformation_,
                    unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction> mEigenfunction_)
            : transformation(transformation_), mEigenfunction(std::move(mEigenfunction_)) {
    }

    Eigen::Array<Scalar, 2, 1> operator()(const Scalar &r) const override {
        Y<Scalar> y;
        Scalar x = transformation->r2x(r);
        y.y() = (*mEigenfunction)(x);
        return transformation->y2z(x, y).y();
    };

    Eigen::Array<Scalar, Eigen::Dynamic, 2>
    operator()(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &r) const override {
        Eigen::Array<Scalar, Eigen::Dynamic, 1> x = r.unaryExpr(
                [&](const Scalar &r_) {
                    return transformation->r2x(r_);
                });

        auto z = (*mEigenfunction)(x);
        Y<Scalar> y;
        for (Eigen::Index i = 0; i < r.rows(); ++i) {
            y.y() = z.row(i);
            z.row(i) = transformation->y2z(x(i), y).y();
        }
        return z;
    };
};

template<typename Scalar>
vector<tuple<int, Scalar, unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction>>>
SturmLiouville<Scalar>::eigenpairsByIndex(
        int Imin, int Imax, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right) const {
    auto schrodingerEigenpairs = matslise.eigenpairsByIndex(
            Imin, Imax,
            transformation.z2y(transformation.rDomain().min(), left),
            transformation.z2y(transformation.rDomain().max(), right));


    vector<tuple<int, Scalar, unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction>>>
            eigenpairs;
    eigenpairs.reserve(schrodingerEigenpairs.size());

    for (auto &iEf: schrodingerEigenpairs) {

        eigenpairs.emplace_back(
                get<0>(iEf),
                get<1>(iEf),
                std::make_unique<SLEigenfunction<Scalar>>(&transformation, std::move(get<2>(iEf)))
        );
    }

    return eigenpairs;
};

#include "./instantiate.h"
