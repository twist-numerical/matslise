#ifndef MATSLISE_LOBATTO_H
#define MATSLISE_LOBATTO_H

#include "eigen.h"

using namespace Eigen;

namespace lobatto {
    ArrayXd grid(const ArrayXd &x);

    double quadrature(const ArrayXd &x, const ArrayXd &f);


    template<int n>
    class Dim {
    };

    template<>
    class Dim<1> {
    public:
        typedef ArrayXd grid;
        typedef std::function<double(double)> function;
    };

    template<>
    class Dim<2> {
    public:
        typedef ArrayXXd grid;
        typedef std::function<double(double, double)> function;
    };


    template<int n>
    typename Dim<n>::grid apply(const ArrayXd grid[n], const typename Dim<n>::function &f);

    template<int n>
    double multi_quadrature(const ArrayXd x[n], const typename Dim<n>::grid &f);
}

#endif