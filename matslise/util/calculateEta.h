#ifndef SCHRODINGER_CALCULATE_ETA_H
#define SCHRODINGER_CALCULATE_ETA_H

#include "eigen.h"

template<typename Scalar>
struct CalculateEtaData {
    static const int taylor_degree;
    static const Scalar taylor_eta8[];
    static const Scalar taylor_eta9[];
};

template<>
inline const int CalculateEtaData<double>::taylor_degree = 7;

#ifdef MATSLISE_long_double
template<>
inline const int CalculateEtaData<long double>::taylor_degree = 11;
#endif


#ifdef MATSLISE_float128
#include <boost/multiprecision/float128.hpp>

template<>
inline const int CalculateEtaData<boost::multiprecision::float128>::taylor_degree = 15;
#endif

template<typename Scalar, int n>
Eigen::Array<Scalar, n, 1> calculateEta(Scalar Z) {
    Eigen::Array<Scalar, n, 1> eta;
    if (abs(Z) < 0.5) {
        Scalar e9 = 0, e8 = 0;
        Scalar z = 1;
        for (int i = 0; i <= CalculateEtaData<Scalar>::taylor_degree; ++i, z *= Z) {
            e8 += z * CalculateEtaData<Scalar>::taylor_eta8[i];
            e9 += z * CalculateEtaData<Scalar>::taylor_eta9[i];
        }

        if constexpr(n > 10)
            eta[10] = e9;
        if constexpr(n > 9)
            eta[9] = e8;

        Scalar eta_prev = e9;
        Scalar eta_curr = e8;
        for (int i = 8; i >= 0; --i) {
            Scalar tmp = eta_curr;
            eta_curr = Z * eta_prev + (2 * i + 1) * eta_curr;
            eta_prev = tmp;
            if (i < n)
                eta[i] = eta_curr;
        }

    } else {
        if (Z > 0) {
            if (Z > 500) {
                // The calculated results will probably be inaccurate
                // std::    cerr << ("Z > 500") << std::endl;
            }
            Scalar sZ = sqrt(Z);
            eta[0] = cosh(sZ);
            eta[1] = sinh(sZ) / sZ;
        } else {
            Scalar sZ = sqrt(-Z);
            eta[0] = cos(sZ);
            eta[1] = sin(sZ) / sZ;
        }

        for (int i = 2; i < n; ++i) {
            eta[i] = (eta[i - 2] - (2 * i - 3) * eta[i - 1]) / Z;
        }
    }

    return eta;
}

template<typename Scalar, int n>
Eigen::Array<Scalar, n, Eigen::Dynamic> calculateEta(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &Z) {
    Eigen::Array<Scalar, n, Eigen::Dynamic> eta(n, Z.size());

    for (int i = 0; i < Z.size(); ++i)
        eta.col(i) = calculateEta<Scalar, n>(Z[i]);

    return eta;
}

#endif