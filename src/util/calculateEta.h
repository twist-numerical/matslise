#ifndef SCHRODINGER_CALCULATE_ETA_H
#define SCHRODINGER_CALCULATE_ETA_H

#include "eigen.h"

template<typename Scalar>
struct CalculateEtaData {
    static const int taylor_degree;
    static const Scalar taylor_eta8[];
    static const Scalar taylor_eta9[];
};


template<typename Scalar>
Scalar *calculateEta(Scalar Z, int etaCount) {
    Scalar *eta;
    if (abs(Z) < 0.5) {
        Scalar e9 = 0, e8 = 0;
        Scalar z = 1;
        for (int i = 0; i <= CalculateEtaData<Scalar>::taylor_degree; ++i, z *= Z) {
            e8 += z * CalculateEtaData<Scalar>::taylor_eta8[i];
            e9 += z * CalculateEtaData<Scalar>::taylor_eta9[i];
        }

        eta = new Scalar[11]{0, 0, 0, 0, 0, 0, 0, 0, 0, e8, e9};
        for (int i = 8; i >= 0; --i)
            eta[i] = Z * eta[i + 2] + (2 * i + 1) * eta[i + 1];

    } else {
        eta = new Scalar[etaCount];

        if (Z > 0) {
            if (Z > 500) {
                // The calculated results will probably be inaccurate
                // std::cerr << ("Z > 500") << std::endl;
            }
            Scalar sZ = sqrt(Z);
            eta[0] = cosh(sZ);
            eta[1] = sinh(sZ) / sZ;
        } else {
            Scalar sZ = sqrt(-Z);
            eta[0] = cos(sZ);
            eta[1] = sin(sZ) / sZ;
        }

        for (int i = 2; i < etaCount; ++i) {
            eta[i] = (eta[i - 2] - (2 * i - 3) * eta[i - 1]) / Z;
        }
    }

    return eta;
}

template<typename Scalar>
Eigen::Array<Scalar, Eigen::Dynamic, 1> *calculateEta(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &Z, int etaCount) {
    auto *eta = new Eigen::Array<Scalar, Eigen::Dynamic, 1>[etaCount];
    for (int j = 0; j < etaCount; ++j)
        eta[j].resize(Z.size(), 1);

    for (int i = 0; i < Z.size(); ++i) {
        Scalar *eta_i = calculateEta<Scalar>(Z[i], etaCount);

        for (int j = 0; j < etaCount; ++j)
            eta[j][i] = eta_i[j];

        delete[] eta_i;
    }

    return eta;
}

#endif