//
// Created by toon on 5/19/20.
//

#ifndef MATSLISE_SECTOR_CORRECTION_POTENTIAL_H
#define MATSLISE_SECTOR_CORRECTION_POTENTIAL_H

#endif //MATSLISE_SECTOR_CORRECTION_POTENTIAL_H

template<typename Scalar>
std::complex<Scalar> sinhc(const std::complex<Scalar> &x) {
    if (abs(x) < .5) {
        // https://en.wikipedia.org/wiki/Sinhc_function#Pad%C3%A9_approximation
        std::complex<Scalar> x2 = x * x;
        return (Scalar(1) + x2 * (Scalar(53272705) / Scalar(360869676) +
                                  x2 * (Scalar(38518909) / Scalar(7217393520) +
                                        x2 * (Scalar(269197963) / Scalar(3940696861920) +
                                              x2 * Scalar(4585922449) / Scalar(15605159573203200))))
               ) / (Scalar(1) + x2 * (Scalar(-2290747) / Scalar(120289892) +
                                      x2 * (Scalar(1281433) / Scalar(7217393520) +
                                            x2 * (Scalar(-560401) / Scalar(562956694560) +
                                                  x2 * Scalar(1029037) / Scalar(346781323848960)))));
    } else
        return sinh(x) / x;
}

template<typename Scalar>
Eigen::Array<std::complex<Scalar>, MATSLISE_HMAX_delta * 2 - 1, 2>
eta_integrals(const Scalar &delta, const std::complex<Scalar> theta) {
    Eigen::Array<std::complex<Scalar>, MATSLISE_HMAX_delta * 2 - 1, 2> integrals = Eigen::Array<Scalar,
            MATSLISE_HMAX_delta * 2 - 1, 2>::Zero();

    std::cout << "θδ(" << theta << ", " << delta << ")" << std::endl;
    if (abs(theta) < Scalar(5) || true) {
        Scalar d = std::pow(delta, MATSLISE_HMAX_delta * 2 - 3);
        for (int i = MATSLISE_HMAX_delta * 2 - 4; i >= 0; --i) {
            // d = delta**(i+1)
            integrals(i, 0) = (d * cosh(delta * theta) - theta * theta * integrals(i + 2, 1)) / Scalar(i + 1);
            if (i > 0)
                integrals(i, 1) = (d * sinhc(delta * theta) - integrals(i, 0)) / Scalar(i);
            d /= delta;
        }
    } else {
        integrals(0, 0) = delta * sinhc(delta * theta);
        integrals(1, 1) = (cosh(theta * delta) - Scalar(1)) / (theta * theta);

        Scalar d = 1;
        for (int i = 1; i < 2 * MATSLISE_HMAX_delta - 1; ++i) {
            if (i > 1)
                integrals(i, 1) = (d * cosh(delta * theta) - Scalar(i - 1) * integrals(i - 2, 0)) / (theta * theta);
            d *= delta;
            integrals(i, 0) = d * delta * sinhc(delta * theta) - Scalar(i) * integrals(i, 1);
        }
    }
    return integrals;
}


template<typename Scalar>
Eigen::Array<Scalar, MATSLISE_N, 1> vbar_formulas(
        Eigen::Array<Scalar, MATSLISE_HMAX_delta, 2> y1,
        Eigen::Array<Scalar, MATSLISE_HMAX_delta, 2> y2,
        const Scalar delta, const Scalar dZ1, const Scalar dZ2) {
    typedef Eigen::Matrix<Scalar, MATSLISE_HMAX_delta, MATSLISE_HMAX_delta> MatrixHs;
    typedef Eigen::Array<std::complex<Scalar>, (MATSLISE_HMAX_delta - 1) * 2 + 1, 1> Array2Hs;

    Eigen::Array<std::complex<Scalar>, MATSLISE_N, 1> quadratures = Eigen::Array<Scalar, MATSLISE_N, 1>::Zero();

    std::complex<Scalar> theta1 = sqrt(std::complex(dZ1));

    for (auto theta2 : std::array<std::complex<Scalar>, 2>{-sqrt(std::complex(dZ2)), sqrt(std::complex(dZ2))}) {
        std::complex<Scalar> theta = theta1 + theta2;
        Array2Hs xi = Array2Hs::Zero(); //eta_-1
        Array2Hs eta = Array2Hs::Zero(); //eta_0
        for (int i = 0; i < MATSLISE_HMAX_delta; ++i) {
            for (int j = 0; j < MATSLISE_HMAX_delta; ++j) {
                xi[i + j] += y1(i, 0) * y2(j, 0) / Scalar(2);
                if (i > 1 && j > 1) {
                    xi[i + j - 2] += y1(i, 1) * y2(j, 1) / Scalar(2) / (theta1 * theta2);
                }

                if (j > 0)
                    eta[i + j] += y1(i, 0) * y2(j, 1) / Scalar(2) * (theta / theta2);
                if (i > 0)
                    eta[i + j] += y1(i, 1) * y2(j, 0) / Scalar(2) * (theta / theta1);
            }
        }
        Eigen::Array<std::complex<Scalar>, 2 * MATSLISE_HMAX_delta - 1, 2> integrals = eta_integrals(delta, theta);

        for (int i = 0; i < MATSLISE_N; ++i) {
            for (int j = 0; i + j < 2 * MATSLISE_HMAX_delta - 1; ++j) {
                quadratures[i] += (integrals(i + j, 0) * xi[j]);
                if (j > 0)
                    quadratures[i] += (integrals(i + j, 1) * eta[j]);
            }
        }

        std::cout << "\n --- y1, y2 " << theta * delta << std::endl;
        std::cout << "\n" << y1.transpose() << std::endl;
        std::cout << "\n" << y2.transpose() << std::endl;
        std::cout << "\n --- xi, eta " << theta * delta << std::endl;
        std::cout << "\n" << integrals.transpose() << std::endl;
        std::cout << "\n" << xi.transpose() << std::endl;
        std::cout << "\n" << eta.transpose() << std::endl;
    }

    std::cout << "\n --- quad " << std::endl;
    std::cout << "\n" << quadratures.transpose() << std::endl;

    return quadratures.real();
}