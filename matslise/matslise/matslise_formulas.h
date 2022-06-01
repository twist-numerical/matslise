#include "../matslise.h"

template<typename Scalar>
void calculate_tcoeff_matrix(
        const Scalar &h,
        const std::array<Scalar, MATSLISE_N> &vs,
        Eigen::Array<Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta, Eigen::DontAlign> &tDelta,
        Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign> *tH);
