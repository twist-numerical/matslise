#ifndef MATSLISE_TEST_PROBLEM_H
#define MATSLISE_TEST_PROBLEM_H

#include "../../matslise/matslise.h"
#include "../test.h"


template<typename Scalar>
Eigen::Array<Scalar, Eigen::Dynamic, 1> lobatto_grid(Eigen::Index intervals, Scalar min, Scalar max) {
    Eigen::Array<Scalar, Eigen::Dynamic, 1> grid(3 * intervals + 1);
    Scalar h = (max - min) / intervals;
    Scalar x1 = Scalar{.5} * h * (Scalar{1} - Scalar{1} / sqrt(Scalar{5}));
    Scalar x2 = h - x1;
    Scalar a = min;
    for (Eigen::Index i = 0; i < intervals; ++i) {
        grid.middleRows(3 * i, 3) << a, a + x1, a + x2;
        a += h;
    }
    grid[3 * intervals] = max;
    return grid;
}

template<typename Scalar>
Scalar lobatto_quadrature(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x,
                          const Eigen::Array<Scalar, Eigen::Dynamic, 1> &f) {
    Eigen::Index n = x.size();
    Scalar result = 0;

    for (Eigen::Index i = 0; i < n - 1; i += 3)
        result += (x[i + 3] - x[i]) / 2 * ((f[i] + f[i + 3]) / 6 + (f[i + 1] + f[i + 2]) * 5 / 6);

    return result;
}

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
        REQUIRE_THAT(iE.second, WithinAbs(correct[j], tolerance) || WithinRel(correct[j], tolerance));
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
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xs = lobatto_grid(501, problem.domain.min, problem.domain.max);

    std::vector<Eigen::Array<Scalar, Eigen::Dynamic, 1>> evaluated;
    for (const auto &iE: eigenvalues) {
        INFO("Checking eigenvalue " << iE.first << ": " << iE.second);

        Eigen::Array<Scalar, Eigen::Dynamic, 2> f_xs = (*problem.eigenfunction(
                iE.second, left, right, iE.first))(xs);
        std::unique_ptr<typename matslise::AbstractMatslise<Scalar>::Eigenfunction> f
                = problem.eigenfunction(iE.second, left, right, iE.first);

        for (int i = 0; i < xs.rows(); ++i) {
            CHECK_THAT(f_xs(i, 0), WithinAbs((*f)(xs[i])(0), tolerance));
        }

        evaluated.emplace_back(f_xs.col(0));
    }

    {
        INFO("Checking orthonormality");
        for (const auto &f1: evaluated)
            for (const auto &f2: evaluated)
                CHECK_THAT(lobatto_quadrature<Scalar>(xs, f1 * f2), WithinAbs<Scalar>(&f1 == &f2 ? 1 : 0, tolerance));
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

    for (int offset: (int[]) {0, n / 6, n / 3, n / 2, 2 * n / 3}) {
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
