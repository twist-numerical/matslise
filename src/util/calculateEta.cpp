//
// Created by toon on 6/26/18.
//

#include <iostream>
#include "calculateEta.h"

double *calculateEta(double Z, int etaCount) {
    double *eta;
    if (fabs(Z) < 0.5) {
        static double eta9[] = {1.527349308567059e-009, 0.36365459727787e-10, 0.00395276736172e-10,
                                0.00002635178241e-10, 0.00000012199899e-10, 0.00000000042069e-10,
                                0.00000000000113e-10, 0};
        static double eta8[] = {2.901963686277412e-008, 0.76367465428353e-9, 0.00909136493195e-9,
                                0.00006587945603e-9, 0.00000032939728e-9, 0.00000000121999e-9,
                                0.00000000000351e-9, 0.00000000000001e-9};

        double e9 = 0, e8 = 0;
        double z = 1;
        for (int i = 0; i < 8; ++i, z *= Z) {
            e9 += z * eta9[i];
            e8 += z * eta8[i];
        }

        eta = new double[11]{0, 0, 0, 0, 0, 0, 0, 0, 0, e8, e9};
        for (int i = 8; i >= 0; --i)
            eta[i] = Z * eta[i + 2] + (2 * i + 1) * eta[i + 1];

    } else {
        eta = new double[etaCount];

        if (Z > 0) {
            if (Z > 500) {
                std::cerr << ("Z > 500") << std::endl;
            }
            double sZ = sqrt(Z);
            eta[0] = cosh(sZ);
            eta[1] = sinh(sZ) / sZ;
        } else {
            double sZ = sqrt(-Z);
            eta[0] = cos(sZ);
            eta[1] = sin(sZ) / sZ;
        }

        for (int i = 2; i < etaCount; ++i) {
            eta[i] = (eta[i - 2] - (2 * i - 3) * eta[i - 1]) / Z;
        }
    }

    return eta;
}

MatrixXd *calculateEta(const VectorXd &Z, int n, int etaCount) {
    VectorXd *eta = new VectorXd[etaCount];
    for (int j = 0; j < etaCount; ++j)
        eta[j] = VectorXd::Zero(n);

    for (int i = 0; i < Z.size(); ++i) {
        double *eta_i = calculateEta(Z[i], etaCount);

        for (int j = 0; j < etaCount; ++j)
            eta[j][i] = eta_i[j];

        delete[] eta_i;
    }

    MatrixXd *eta_mat = new MatrixXd[etaCount];
    for (int j = 0; j < etaCount; ++j)
        eta_mat[j] = eta[j].asDiagonal();
	delete[] eta;
    return eta_mat;
}