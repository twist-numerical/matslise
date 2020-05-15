//
// Created by toon on 4/20/20.
//

#ifndef MATSLISE_TEST_PROBLEM_H
#define MATSLISE_TEST_PROBLEM_H

#include "../../src/matslise.h"
#include "../../src/util/lobatto.h"
#include "../catch.hpp"

template<typename Scalar>
void testEigenvaluesByIndex(
        const matslise::AbstractMatslise<Scalar> &problem,
        const matslise::Y<Scalar> left,
        const matslise::Y<Scalar> right,
        const std::vector<Scalar> &correct,
        const Scalar tolerance,
        const int offset,
        const int count) {


    auto eigenvalues = problem.eigenvaluesByIndex(offset, offset + count, left, right);
    int j = offset;
    for (const auto &iE: eigenvalues) {
        REQUIRE(iE.first == j);
        REQUIRE(Approx(iE.second).margin(tolerance) == correct[j]);
        ++j;
    }
}

template<typename Scalar>
void testOrthogonality(
        const matslise::AbstractMatslise<Scalar> &problem,
        const matslise::Y<Scalar> left,
        const matslise::Y<Scalar> right,
        const std::vector<std::pair<int, Scalar>> &eigenvalues,
        const Scalar &tolerance) {
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xs = lobatto::grid<Scalar>(
            Eigen::Array<Scalar, Eigen::Dynamic, 1>::LinSpaced(
                    501, problem.domain.min, problem.domain.max));

    std::vector<Eigen::Array<Scalar, Eigen::Dynamic, 1>> evaluated;
    for (const auto &iE:  eigenvalues) {
        INFO("Checking eigenvalue " << iE.first << ": " << iE.second);

        Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1> f_xs = problem.eigenfunction(
                iE.second, left, right, iE.first)(xs);
        std::function<matslise::Y<Scalar>(Scalar)> f = problem.eigenfunction(iE.second, left, right, iE.first);

        for (int i = 0; i < xs.rows(); ++i) {
            CHECK(Approx(f_xs[i].y[0]).margin(tolerance) == f(xs[i]).y[0]);
        }

        evaluated.emplace_back(f_xs.template unaryExpr<std::function<Scalar(const matslise::Y<Scalar> &)>>(
                [&](const matslise::Y<Scalar> &y) { return y.y[0]; }));
    }

    {
        INFO("Checking orthonormality")
        for (const auto &f1 : evaluated)
            for (const auto &f2 : evaluated)
                CHECK(Approx(lobatto::quadrature<Scalar>(xs, f1 * f2)).margin(tolerance) == (&f1 == &f2 ? 1 : 0));
    }
}

template<typename Scalar>
void testProblem(
        const matslise::AbstractMatslise<Scalar> &problem,
        const matslise::Y<Scalar> left,
        const matslise::Y<Scalar> right,
        const std::vector<Scalar> &correct,
        const Scalar &tolerance, const Scalar &toleranceEigenfunctions) {

    int n = (int) correct.size();

    for (int offset : (int[]) {0, n / 6, n / 3, n / 2, 2 * n / 3}) {
        testEigenvaluesByIndex(problem, left, right, correct, tolerance, offset, n - offset);
    }

    testOrthogonality(problem, left, right, problem.eigenvaluesByIndex(0, correct.size(), left, right),
                      toleranceEigenfunctions);

}

template<typename Scalar>
void testProblem(
        const matslise::AbstractMatslise<Scalar> &problem,
        const matslise::Y<Scalar> left,
        const matslise::Y<Scalar> right,
        const std::vector<Scalar> &correct,
        const Scalar &tolerance) {
    testProblem(problem, left, right, correct, tolerance, tolerance);
}

#endif //MATSLISE_TEST_PROBLEM_H
