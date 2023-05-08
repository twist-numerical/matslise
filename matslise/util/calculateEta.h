#ifndef SCHRODINGER_CALCULATE_ETA_H
#define SCHRODINGER_CALCULATE_ETA_H

#include "eigen.h"
#include <math.h>

template<typename Scalar, int n>
Eigen::Array<Scalar, n, 1> calculateEta(Scalar Z) {
    Eigen::Array<Scalar, n, 1> eta;
    if (abs(Z) < 100) {
        int last_eta = (std::abs(int(Z)) + 9) / 4;
        if (last_eta > 150) last_eta = 150;
        if (last_eta < n + 2) last_eta = n + 2;

        Scalar w0 = 1;
        for (int j = 0; j < last_eta; ++j)
            w0 *= 2 * j + 1;
        w0 = 1 / w0;
        // eta_current
        Scalar eta_current = w0;
        {
            Scalar f = w0;
            for (int q2 = 2; q2 < 60; q2 += 2) {
                f *= Z / (q2 * (q2 + 2 * last_eta - 1));
                eta_current += f;
            }
        }

        // eta_next
        w0 /= (2 * last_eta + 1);
        Scalar eta_next = w0;
        {
            Scalar f = w0;
            for (int q2 = 2; q2 < 60; q2 += 2) {
                f *= Z / (q2 * (q2 + 2 * last_eta + 1));
                eta_next += f;
            }
        }

        // eta_(last_eta) and eta_(last_eta-1)
        for (int i = last_eta - 1; i >= 0; --i) {
            Scalar tmp = eta_current;
            eta_current = (2 * i + 1) * eta_current + Z * eta_next;
            eta_next = tmp;
            if (i < n)
                eta[i] = eta_current;
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
