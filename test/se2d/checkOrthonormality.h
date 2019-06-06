#ifndef MATSLISE_CHECKORTHONORMALITY_H
#define MATSLISE_CHECKORTHONORMALITY_H

#include "../catch.hpp"
#include "../../src/matslise.h"
#include "../../src/util/lobatto.h"

using namespace std;
using namespace Eigen;
using namespace matslise;

template<typename doubleIterator>
void checkOrthonormality(const SEnD<2> &p, const doubleIterator &begin, const doubleIterator &end) {
    int n = 71;
    ArrayXd xy[2] = {
            lobatto::grid(ArrayXd::LinSpaced(n, p.domain.getMin(0), p.domain.getMax(0))),
            lobatto::grid(ArrayXd::LinSpaced(n, p.domain.getMin(1), p.domain.getMax(1)))
    };

    vector<ArrayXXd> eigenfunctions;
    for (auto i = begin; i != end; ++i) {
        vector<ArrayXXd> *fs = p.computeEigenfunction(*i, xy);
        for (const ArrayXXd &f : *fs)
            eigenfunctions.push_back(f);
        delete fs;
    }

    for (auto i = eigenfunctions.begin(); i != eigenfunctions.end(); ++i)
        for (auto j = eigenfunctions.begin(); j != eigenfunctions.end(); ++j) {
            CHECKED_ELSE(Approx(lobatto::multi_quadrature<2>(xy, *i * *j)).margin(1e-1) == (i == j ? 1 : 0)) {
                auto l = begin;
                for (auto k = eigenfunctions.begin(); k != eigenfunctions.end(); ++k) {
                    if (k == i)
                        break;
                    ++l;
                }
                FAIL(*l);
            }
        }

}


#endif //MATSLISE_CHECKORTHONORMALITY_H
