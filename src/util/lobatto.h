#include "../Matrix2D.h"

using namespace Eigen;

namespace lobatto {
    ArrayXd grid(const ArrayXd &x);

    double quadrature(const ArrayXd &x, const ArrayXd &f);
}