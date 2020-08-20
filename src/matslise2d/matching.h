#ifndef MATSLISE_MATCHING_H
#define MATSLISE_MATCHING_H

#include "../util/eigen.h"
#include "../matslise.h"
#include <tuple>
#include <map>

using namespace matslise;
using namespace Eigen;

template<typename Scalar>
std::pair<Matrix<Scalar, Dynamic, Dynamic>, Matrix<Scalar, Dynamic, Dynamic>> errorMatrix(
        const Y<Scalar, Dynamic> &left, const Y<Scalar, Dynamic> &right) {
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXs;
    ColPivHouseholderQR<MatrixXs> left_solver(left.getY(0).transpose());
    ColPivHouseholderQR<MatrixXs> right_solver(right.getY(0).transpose());
    MatrixXs Ul = left_solver.solve(left.getY(1).transpose()).transpose();
    MatrixXs Ur = right_solver.solve(right.getY(1).transpose()).transpose();
    return {
            Ul - Ur,
            left_solver.solve((left.getdY(1) - Ul * left.getdY(0)).transpose()).transpose()
            - right_solver.solve((right.getdY(1) - Ur * right.getdY(0)).transpose()).transpose()
    };
}

template<typename Scalar, bool withVectors>
std::vector<typename std::conditional<withVectors,
        std::tuple<Scalar, Scalar, Matrix<Scalar, Dynamic, 1>>,
        std::pair<Scalar, Scalar>>::type>
eigenvaluesWithDerivatives(const Matrix<Scalar, Dynamic, Dynamic> &error,
                           const Matrix<Scalar, Dynamic, Dynamic> &derivative) {
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXs;
    typedef Matrix<Scalar, Dynamic, 1> VectorXs;
    assert(error.rows() == error.cols() && error.rows() == derivative.rows() && error.cols() == derivative.cols());
    Index N = error.rows();

    EigenSolver<MatrixXs> solver(N);

    solver.compute(error, true);

    std::multimap<Scalar, int> rightMap;
    Array<Scalar, Dynamic, 1> eigenvaluesRight = solver.eigenvalues().array().real();
    for (int i = 0; i < N; ++i)
        rightMap.insert({eigenvaluesRight[i], i});
    // We assume that .eigenvalues() and .pseudoEigenvectors() return elements in the same order
    // We assume that .pseudoEigenvectors() will fix some issues with complex eigenvectors.
    MatrixXs right = solver.pseudoEigenvectors();

    std::multimap<Scalar, int> leftMap;
    solver.compute(error.transpose(), true);
    Array<Scalar, Dynamic, 1> eigenvaluesLeft = solver.eigenvalues().array().real();
    for (int i = 0; i < N; ++i)
        leftMap.insert({eigenvaluesLeft[i], i});
    MatrixXs left = solver.pseudoEigenvectors().transpose();

    std::vector<typename std::conditional<withVectors,
            std::tuple<Scalar, Scalar, VectorXs>,
            std::pair<Scalar, Scalar>>::type> errors;
    for (auto leftI = leftMap.begin(), rightI = rightMap.begin();
         leftI != leftMap.end() && rightI != rightMap.end();
         ++leftI, ++rightI) {
        int &li = leftI->second;
        int &ri = rightI->second;
        Scalar e = eigenvaluesLeft[li];
        Scalar de = (left.row(li) * derivative * right.col(ri) /
                     (left.row(li) * right.col(ri)))[0];
        if constexpr (withVectors) {
            errors.emplace_back(e, de, right.col(ri));
        } else {
            errors.emplace_back(e, de);
        }
    }

    return errors;
}

template<typename Scalar>
Matrix<Scalar, Dynamic, Dynamic>
getKernel(const Y<Scalar, Dynamic> &left, const Y<Scalar, Dynamic> &right, const Scalar &tolerance) {
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXs;
    typedef Matrix<Scalar, Dynamic, 1> VectorXs;
    std::pair<MatrixXs, MatrixXs> error_matrix = errorMatrix<Scalar>(left, right);
    std::vector<VectorXs> kernel_list;
    for (auto &error : eigenvaluesWithDerivatives<Scalar, true>(error_matrix.first, error_matrix.second)) {
        if (abs(std::get<0>(error) / std::get<1>(error)) < tolerance) {
            kernel_list.push_back(std::move(std::get<2>(error)));
        }
    }
    MatrixXs kernel(left.getN(), kernel_list.size());
    for (Index i = 0; i < static_cast<Index>(kernel_list.size()); ++i)
        kernel.col(i) = kernel_list[i];
    return kernel;
}


#endif //MATSLISE_MATCHING_H
