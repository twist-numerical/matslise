#include "matslise_formulas.h"

void calculate_tcoeff_matrix(double h, double *vs, Array2D<Matrix2D<double>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta> &tDelta, Matrix2D<double> *tH) {
    double &v1 = vs[1];
    double &v2 = vs[2];
    double &v3 = vs[3];
    double &v4 = vs[4];
    double &v5 = vs[5];
    double &v6 = vs[6];
    tDelta[0][0] = {(1), 0, 0, (1)}; 
    tDelta[0][1] = {0, 0, ((((-0.5*v1))+(((0.5*v2))+(((-0.5*v3))+((0.5*v4))*h)*h)*h)*h), 0}; 
    tDelta[0][2] = {0, 0, (((0.5*v1))+(((-1.5*v2))+(((3.0*v3))+((-5.0*v4))*h)*h)*h), 0}; 
    tDelta[0][3] = {0, 0, ((v2)+(((-5.0*v3))+((15.0*v4))*h)*h), 0}; 
    tDelta[0][4] = {0, 0, (((2.5*v3))+((-17.5*v4))*h), 0}; 
    tDelta[0][5] = {0, 0, ((7.0*v4)), 0}; 
    tDelta[0][6] = {0, 0, 0, 0}; 
    tH[0] = {(1), 0, 0, (1)}; 
    tDelta[1][0] = {0, 0, 0, 0}; 
    tDelta[1][1] = {0, (1), ((((-0.5*v1))+(((0.5*v2))+(((-0.5*v3))+((0.5*v4))*h)*h)*h)*h), 0}; 
    tDelta[1][2] = {((((-0.5*v1))+(((0.5*v2))+(((-0.5*v3))+((0.5*v4))*h)*h)*h)*h), 0, (((0.5*v1))+(((-1.5*v2))+(((3.0*v3))+((-5.0*v4))*h)*h)*h), ((((-0.5*v1))+(((0.5*v2))+(((-0.5*v3))+((0.5*v4))*h)*h)*h)*h)}; 
    tDelta[1][3] = {(((0.5*v1))+(((-1.5*v2))+(((3.0*v3))+((-5.0*v4))*h)*h)*h), 0, (((1.5*v2))+(((-7.5*v3))+(((0.125*(v1*v1))+(22.5*v4))+((-0.25*(v1*v2)))*h)*h)*h), (((0.5*v1))+(((-1.5*v2))+(((3.0*v3))+((-5.0*v4))*h)*h)*h)}; 
    tDelta[1][4] = {((v2)+(((-5.0*v3))+((15.0*v4))*h)*h), 0, (((5.0*v3))+(((-0.25*(v1*v1))+(-35.0*v4))+((v1*v2))*h)*h), ((v2)+(((-5.0*v3))+((15.0*v4))*h)*h)}; 
    tDelta[1][5] = {(((2.5*v3))+((-17.5*v4))*h), 0, (((0.125*(v1*v1))+(17.5*v4))+((-1.25*(v1*v2)))*h), (((2.5*v3))+((-17.5*v4))*h)}; 
    tDelta[1][6] = {((7.0*v4)), 0, ((0.5*(v1*v2))), ((7.0*v4))}; 
    tH[1] = {0, ((1)*h), ((((((0.5*v2))+((((0.5*v4))+(((0.5*v6))*h)*h)*h)*h)*h)*h)*h), 0}; 
    tDelta[2][0] = {0, 0, 0, 0}; 
    tDelta[2][1] = {0, 0, 0, 0}; 
    tDelta[2][2] = {0, 0, 0, 0}; 
    tDelta[2][3] = {(((-0.5*v1))+(((1.5*v2))+(((-3.0*v3))+((5.0*v4))*h)*h)*h), ((((-0.5*v1))+(((0.5*v2))+((-0.5*v3))*h)*h)*h), (((-1.5*v2))+(((7.5*v3))+(((0.125*(v1*v1))+(-22.5*v4))+((-0.25*(v1*v2)))*h)*h)*h), (((0.5*v1))+(((-1.5*v2))+(((3.0*v3))+((-5.0*v4))*h)*h)*h)}; 
    tDelta[2][4] = {(((-1.5*v2))+(((7.5*v3))+((0.125*(v1*v1))+(-22.5*v4))*h)*h), (((0.5*v1))+(((-1.5*v2))+((3.0*v3))*h)*h), (((-7.5*v3))+(((-0.25*(v1*v1))+(52.5*v4))+((v1*v2))*h)*h), (((1.5*v2))+(((-7.5*v3))+((0.125*(v1*v1))+(22.5*v4))*h)*h)}; 
    tDelta[2][5] = {(((-5.0*v3))+((-0.25*(v1*v1))+(35.0*v4))*h), ((v2)+((-5.0*v3))*h), (((0.08333333333333333*(v1*v1))+(-35.0*v4))+((-1.25*(v1*v2)))*h), (((5.0*v3))+((-0.25*(v1*v1))+(-35.0*v4))*h)}; 
    tDelta[2][6] = {((0.125*(v1*v1))+(-17.5*v4)), ((2.5*v3)), ((0.5*(v1*v2))), ((0.125*(v1*v1))+(17.5*v4))}; 
    tH[2] = {((((((-0.5*v1))+((((-0.5*v3))+(((-0.5*v5))*h)*h)*h)*h)*h)*h)*h), 0, ((((((-1.5*v2))+((((-5.0*v4)+(-0.041666666666666664*(v1*v1)))+((((-10.5*v6)+(-0.025*(v2*v2)))+((((-0.017857142857142856*(v3*v3)))+((((-0.013888888888888888*(v4*v4)))+(((-0.011363636363636364*(v5*v5)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h), ((((((0.5*v1))+((((0.5*v3))+(((0.5*v5))*h)*h)*h)*h)*h)*h)*h)}; 
    tH[3] = {((((((((2.5*v3))+(((-0.041666666666666664*(v1*v1)))+(((7.0*v5))+(((-0.025*(v2*v2)))+((((-0.017857142857142856*(v3*v3)))+((((-0.013888888888888888*(v4*v4)))+(((-0.011363636363636364*(v5*v5)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h), ((((((((-0.5*v2))+((((-0.5*v4))+(((-0.5*v6))*h)*h)*h)*h)*h)*h)*h)*h)*h), ((((((((17.5*v4)+(-0.2916666666666667*(v1*v1)))+((((-0.5*(v1*v3))+(94.5*v6)+(-0.15*(v2*v2)))+((((-0.5*(v1*v5))+(-0.25*(v2*v4))+(-0.26785714285714285*(v3*v3)))+((((-0.5*(v3*v5))+(-0.1388888888888889*(v4*v4))+(-0.25*(v2*v6)))+(((-0.26136363636363635*(v5*v5))+(-0.25*(v4*v6)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h), ((((((((-2.5*v3))+(((-0.041666666666666664*(v1*v1)))+(((-7.0*v5))+(((-0.025*(v2*v2)))+((((-0.017857142857142856*(v3*v3)))+((((-0.013888888888888888*(v4*v4)))+(((-0.011363636363636364*(v5*v5)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)}; 
    tH[4] = {((((((((((-31.5*v5)+(0.5*(v1*v2)))+(((0.25*(v1*v3))+(-0.075*(v2*v2)))+(((0.5*(v1*v4))+(0.5*(v2*v3)))+(((0.25*(v1*v5))+(0.07142857142857142*(v3*v3)))+(((0.5*(v2*v5))+(0.5*(v1*v6))+(0.5*(v3*v4)))+(((0.25*(v3*v5))+(-0.041666666666666664*(v4*v4)))+(((0.5*(v4*v5))+(0.5*(v3*v6)))+(((0.09090909090909091*(v5*v5)))+((0.5*(v5*v6)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h), ((((((((((3.5*v4)+(-0.041666666666666664*(v1*v1)))+((((9.0*v6)+(-0.025*(v2*v2)))+((((-0.017857142857142856*(v3*v3)))+((((-0.013888888888888888*(v4*v4)))+(((-0.011363636363636364*(v5*v5)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h), ((((((((((4.0*(v1*v3))+(-346.5*v6)+(1.425*(v2*v2)))+((((10.75*(v1*v5))+(6.5*(v2*v4))+(0.004166666666666667*(v1*v1*v2))+(4.446428571428571*(v3*v3)))+((((15.75*(v3*v5))+(5.833333333333333*(v4*v4))+(12.0*(v2*v6))+(-0.020833333333333332*(v1*v1*v4))+(-0.008928571428571428*(v2*v2*v2))+(0.03214285714285714*(v1*v2*v3)))+(((-0.020833333333333332*(v1*v1*v6))+(0.023809523809523808*(v1*v3*v4))+(-0.0017857142857142857*(v2*v2*v4))+(-0.0017857142857142857*(v2*v3*v3))+(12.340909090909092*(v5*v5))+(17.25*(v4*v6)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h), ((((((((((31.5*v5)+(-0.5*(v1*v2)))+(((0.25*(v1*v3))+(-0.075*(v2*v2)))+(((-0.5*(v1*v4))+(-0.5*(v2*v3)))+(((0.25*(v1*v5))+(0.07142857142857142*(v3*v3)))+(((-0.5*(v2*v5))+(-0.5*(v1*v6))+(-0.5*(v3*v4)))+(((0.25*(v3*v5))+(-0.041666666666666664*(v4*v4)))+(((-0.5*(v4*v5))+(-0.5*(v3*v6)))+((0.09090909090909091*(v5*v5)))*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)*h)}; 
}
