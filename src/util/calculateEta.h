#ifndef SCHRODINGER_CALCULATE_ETA_H
#define SCHRODINGER_CALCULATE_ETA_H

#include "eigen.h"
#include "fmath.h"


using namespace Eigen;

template<typename Scalar>
Scalar *calculateEta(Scalar Z, int etaCount) {
    Scalar *eta;
    if (fmath<Scalar>::abs(Z) < 0.5) {
        static Scalar eta9[] = {1.527349308567059e-009, 0.36365459727787e-10, 0.00395276736172e-10,
                                0.00002635178241e-10, 0.00000012199899e-10, 0.00000000042069e-10,
                                0.00000000000113e-10, 0};
        static Scalar eta8[] = {2.901963686277412e-008, 0.76367465428353e-9, 0.00909136493195e-9,
                                0.00006587945603e-9, 0.00000032939728e-9, 0.00000000121999e-9,
                                0.00000000000351e-9, 0.00000000000001e-9};

        Scalar e9 = 0, e8 = 0;
        Scalar z = 1;
        for (int i = 0; i < 8; ++i, z *= Z) {
            e9 += z * eta9[i];
            e8 += z * eta8[i];
        }

        eta = new Scalar[11]{0, 0, 0, 0, 0, 0, 0, 0, 0, e8, e9};
        for (int i = 8; i >= 0; --i)
            eta[i] = Z * eta[i + 2] + (2 * i + 1) * eta[i + 1];

    } else {
        eta = new Scalar[etaCount];

        if (Z > 0) {
            if (Z > 500) {
                std::cerr << ("Z > 500") << std::endl;
            }
            Scalar sZ = fmath<Scalar>::sqrt(Z);
            eta[0] = fmath<Scalar>::cosh(sZ);
            eta[1] = fmath<Scalar>::sinh(sZ) / sZ;
        } else {
            Scalar sZ = fmath<Scalar>::sqrt(-Z);
            eta[0] = fmath<Scalar>::cos(sZ);
            eta[1] = fmath<Scalar>::sin(sZ) / sZ;
        }

        for (int i = 2; i < etaCount; ++i) {
            eta[i] = (eta[i - 2] - (2 * i - 3) * eta[i - 1]) / Z;
        }
    }

    return eta;
}

template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic> *calculateEta(const Matrix<Scalar, Dynamic, 1> &Z, int n, int etaCount) {
    Matrix<Scalar, Dynamic, 1> *eta = new Matrix<Scalar, Dynamic, 1>[etaCount];
    for (int j = 0; j < etaCount; ++j)
        eta[j] = Matrix<Scalar, Dynamic, 1>::Zero(n);

    for (int i = 0; i < Z.size(); ++i) {
        Scalar *eta_i = calculateEta<Scalar>(Z[i], etaCount);

        for (int j = 0; j < etaCount; ++j)
            eta[j][i] = eta_i[j];

        delete[] eta_i;
    }

    Matrix<Scalar, Dynamic, Dynamic> *eta_mat = new Matrix<Scalar, Dynamic, Dynamic>[etaCount];
    for (int j = 0; j < etaCount; ++j)
        eta_mat[j] = eta[j].asDiagonal();

    delete[] eta;
    return eta_mat;
}

#endif