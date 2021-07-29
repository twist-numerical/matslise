#include <vector>
#include "../matslise.h"

using namespace matslise;
using namespace std;
using namespace Eigen;

template<typename Scalar>
vector<tuple<int, Scalar, typename AbstractMatslise<Scalar>::Eigenfunction>>
AbstractMatslise<Scalar>::eigenpairsByIndex(int Imin, int Imax, const Y<Scalar> &left, const Y<Scalar> &right) const {
    vector<tuple<int, Scalar, AbstractMatslise<Scalar>::Eigenfunction>> result;
    vector<pair<int, Scalar>> eigenvalues = eigenvaluesByIndex(Imin, Imax, left, right);
    result.reserve(eigenvalues.size());
    for (auto &v: eigenvalues) {
        result.emplace_back(v.first, v.second, eigenfunction(v.second, left, right, v.first));
    }
    return result;
}


#include "instantiate.h"
