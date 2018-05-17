//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_LEGENDRE_H
#define SCHRODINGER_LEGENDRE_H


namespace legendre {
    double *getCoefficients(int n, std::function<double(double)> V, double a, double b);
};


#endif //SCHRODINGER_LEGENDRE_H
