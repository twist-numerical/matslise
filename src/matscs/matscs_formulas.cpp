#include "matscs_formulas.h"

void calculate_tcoeff_matrix(int n, double h, MatrixXd *vs, Array2D<Matrix2D<MatrixXd>, MATSCS_ETA, MATSCS_HMAX> &t) {
    MatrixXd zero = MatrixXd::Zero(n, n);
    MatrixXd I = MatrixXd::Identity(n, n);
    MatrixXd &v0 = vs[0];
    MatrixXd &v1 = vs[1];
    MatrixXd &v2 = vs[2];
    t[0][0] = {(I), zero, zero, (I)}; 
    t[0][1] = {zero, zero, ((((-0.5*v1))+((0.5*v2))*h)*h), zero}; 
    t[0][2] = {zero, zero, (((0.5*v1))+((-1.5*v2))*h), zero}; 
    t[0][3] = {zero, zero, (v2), zero}; 
    t[0][4] = {zero, zero, zero, zero}; 
    t[1][0] = {zero, zero, zero, zero}; 
    t[1][1] = {zero, (I), ((((-0.5*v1))+((0.5*v2))*h)*h), zero}; 
    t[1][2] = {((((-0.5*v1))+((0.5*v2))*h)*h), zero, (((0.5*v1))+((-1.5*v2))*h), ((((-0.5*v1))+((0.5*v2))*h)*h)}; 
    t[1][3] = {(((0.5*v1))+((-1.5*v2))*h), zero, (((1.5*v2))+((0.125*(v1*v0))+(-0.125*(v0*v1)))*h), (((0.5*v1))+((-1.5*v2))*h)}; 
    t[1][4] = {(v2), zero, ((-0.08333333333333333*(v1*v0))+(0.08333333333333333*(v0*v1))), (v2)}; 
    t[2][0] = {zero, zero, zero, zero}; 
    t[2][1] = {zero, zero, zero, zero}; 
    t[2][2] = {zero, zero, zero, zero}; 
    t[2][3] = {(((-0.5*v1))+((1.5*v2))*h), (((-0.5*v1))*h), (((-1.5*v2))+((0.125*(v1*v0))+(-0.125*(v0*v1)))*h), (((0.5*v1))+((-1.5*v2))*h)}; 
    t[2][4] = {((-1.5*v2)), ((0.5*v1)), zero, ((1.5*v2))}; 
}
