#ifndef SCHRODINGER_MATSCS_FORMULAS_H
#define SCHRODINGER_MATSCS_FORMULAS_H

#include "../formula_constants.h"
#include "../util/eigen.h"
#include <array>

template<typename Scalar>
void calculate_tcoeff_matrix(
        int n,
        Scalar h,
        const std::array<Eigen::Matrix<Scalar, -1, -1>, MATSCS_N> &vs,
        Eigen::Array<Eigen::Matrix<Scalar, -1, -1>, MATSCS_ETA_delta, MATSCS_HMAX_delta> &tDelta,
        std::array<Eigen::Matrix<Scalar, -1, -1>, MATSCS_ETA_h> &tH);

#endif
