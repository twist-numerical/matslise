#ifndef MATSLISE_FEQUALS_H
#define MATSLISE_FEQUALS_H

#include <random>
#include "../test.h"

template<typename Scalar, typename F1, typename F2>
void fEquals(F1 f1, F2 f2, Scalar xmin, Scalar xmax, Scalar margin = 1e-6) {
    static thread_local std::default_random_engine gen(std::random_device{}());

    Scalar h = (xmax - xmin) / 256;
    std::uniform_real_distribution<Scalar> dis(0, h);
    for (double x = xmin; x + h <= xmax; x += h) {
        double v = x + dis(gen);
        REQUIRE_THAT(f1(v), WithinAbs(f2(v), margin));
    }
}

#endif //MATSLISE_FEQUALS_H
