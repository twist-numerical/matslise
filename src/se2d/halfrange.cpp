#include <iostream>
#include <map>
#include "../matslise.h"
#include "../util/lobatto.h"
#include "../util/find_sector.h"

using namespace Eigen;
using namespace matslise;
using namespace std;


template<typename Scalar>
SE2DHalf<Scalar>::SE2DHalf(const function<Scalar(const Scalar &, const Scalar &)> &V,
                           const matslise::Rectangle<2, Scalar> &domain,
                           const Options2<Scalar> &options) {
    se2d = new SE2D<Scalar>(V, {domain.sub, 0, domain.getMax(1)}, options);
    evenBoundary = Y<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Neumann(se2d->N);
    oddBoundary = Y<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Dirichlet(se2d->N);
}

template<typename Scalar>
SE2DHalf<Scalar>::~SE2DHalf() {
    delete se2d;
}

template<typename Scalar>
void SE2DHalf<Scalar>::setParity(bool even) {
    se2d->y0Left = even ? evenBoundary : oddBoundary;
}

template<typename Scalar>
Scalar SE2DHalf<Scalar>::findEigenvalue(
        const Scalar &guessE, const Scalar &tolerance, int maxIterations, const Scalar &minTolerance) {
    setParity(true);
    Scalar evenE = se2d->findEigenvalue(guessE, tolerance, maxIterations, minTolerance);
    setParity(false);
    Scalar oddE = se2d->findEigenvalue(guessE, tolerance, maxIterations, minTolerance);
    if (abs(evenE - guessE) < abs(oddE - guessE))
        return evenE;
    else
        return oddE;
}

template<typename Scalar>
vector<Scalar> SE2DHalf<Scalar>::findEigenvalues(
        const Scalar &Emin, const Scalar &Emax, const int &initialSteps) {
    setParity(true);
    vector<Scalar> eigenvalues = se2d->findEigenvalues(Emin, Emax, initialSteps);
    setParity(false);
    vector<Scalar> eigenvaluesOdd = se2d->findEigenvalues(Emin, Emax, initialSteps);
    eigenvalues.insert(eigenvalues.end(),
                       make_move_iterator(eigenvaluesOdd.begin()),
                       make_move_iterator(eigenvaluesOdd.end()));
    sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
}

template<typename Scalar>
Scalar SE2DHalf<Scalar>::findFirstEigenvalue() {
    setParity(true);
    return se2d->findFirstEigenvalue();
}

template<typename Scalar>
std::vector<typename SE2DHalf<Scalar>::ArrayXXs>
SE2DHalf<Scalar>::computeEigenfunction(const Scalar &, const ArrayXs &, const ArrayXs &) {
    throw std::logic_error("Function not yet implemented");
}

template<typename Scalar>
std::vector<Scalar> SE2DHalf<Scalar>::computeEigenvaluesByIndex(int, int) {
    throw std::logic_error("Function not yet implemented");
}
// std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

template<typename Scalar>
std::vector<std::function<Scalar(Scalar, Scalar)>> SE2DHalf<Scalar>::eigenfunctionCalculator(const Scalar &) {
    throw std::logic_error("Function not yet implemented");
}


#include "../util/instantiate.h"