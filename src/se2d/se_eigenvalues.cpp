#include "../matslise.h"

using namespace Eigen;
using namespace matslise;
using namespace std;

template<typename Scalar>
vector<Scalar> *SE2D<Scalar>::computeEigenvaluesByIndex(int Imin, int Imax) const {
    // TODO
    (void) Imin;
    (void) Imax;
    auto values = new vector<Scalar>();
    return values;
}

#include "../util/instantiate.h"