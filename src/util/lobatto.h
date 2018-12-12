#ifndef MATSLISE_LOBATTO_H
#define MATSLISE_LOBATTO_H

#include "../Matrix2D.h"

using namespace Eigen;

namespace lobatto {
    ArrayXd grid(const ArrayXd &x);

    double quadrature(const ArrayXd &x, const ArrayXd &f);


    template<int n>
    class Dim {};

    template<>
    class Dim<1> {
    public:
        typedef ArrayXd grid;
    };

    template<>
    class Dim<2> {
    public:
        typedef ArrayXXd grid;
    };

    template<int n>
    double multi_quadrature(const ArrayXd x[n], const typename Dim<n>::grid &f);
}

#endif