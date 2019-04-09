#include "../matscs.h"

void calculate_tcoeff_matrix(int n, double h, MatrixXd *vs, Array2D<MatrixXd, MATSCS_ETA_delta, MATSCS_HMAX_delta> &tDelta, MatrixXd *tH);