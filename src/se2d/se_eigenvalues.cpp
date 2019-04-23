#include "../se2d.h"

using namespace Eigen;
using namespace matslise;
using namespace matslise::SEnD_util;
using namespace std;

template<>
vector<double> *SEnD<2>::computeEigenvaluesByIndex(int Imin, int Imax) const {
    auto values = new vector<double>();

    return values;
}