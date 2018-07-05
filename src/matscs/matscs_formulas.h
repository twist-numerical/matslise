#include "../matslise.h"
#include "../matscs.h"
#include "../Array2D.h"
#include "../Matrix2D.h"

void calculate_tcoeff_matrix(int n, double h, MatrixXd *vs, Array2D<Matrix2D<MatrixXd>, MATSCS_ETA, MATSCS_HMAX> &t);