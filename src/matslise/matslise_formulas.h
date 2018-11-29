#include "../matslise.h"
#include "../Array2D.h"
#include "../Matrix2D.h"

void calculate_tcoeff_matrix(double h, double *vs, Array2D<Matrix2D<double>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &tDelta, Matrix2D<double> *tH);