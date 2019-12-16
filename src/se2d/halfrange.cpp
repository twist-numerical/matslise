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
                           const Options2<Scalar> &options) : domain(domain) {
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
pair<Scalar, Scalar> SE2DHalf<Scalar>::calculateError(
        const Scalar &E, const function<bool(std::pair<Scalar, Scalar>, std::pair<Scalar, Scalar>)> &sorter) {
    setParity(true);
    pair<Scalar, Scalar> errorEven = se2d->calculateError(E, sorter);
    setParity(false);
    pair<Scalar, Scalar> errorOdd = se2d->calculateError(E, sorter);
    return abs(get<0>(errorEven)) <= abs(get<0>(errorOdd)) ? errorEven : errorOdd;
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
vector<typename SE2DHalf<Scalar>::ArrayXXs>
SE2DHalf<Scalar>::computeEigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) {
    Eigen::Index n = y.size();
    for (int i = 1; i < n; ++i)
        if (y[i - 1] > y[i])
            throw runtime_error("SE2DHalf::computeEigenfunction(): y has to be sorted");

    Eigen::Index negatives = 0;
    while(negatives < n && y[negatives] < 0)
        ++negatives;
    cout << negatives << endl;

    ArrayXs yNeg = -y.topRows(negatives).reverse();
    ArrayXs yPos = y.bottomRows(n - negatives);


   // const Scalar SQRT1_2 = sqrt(Scalar(.5));

    vector<ArrayXXs> result;
    for (bool even : {false, true}) {
        setParity(even);
        if (abs(E - se2d->findEigenvalue(E)) < 1e-3) {
            vector<ArrayXXs> fsNeg = se2d->computeEigenfunction(E, x, yNeg);
            vector<ArrayXXs> fsPos = se2d->computeEigenfunction(E, x, yPos);
            for (unsigned int i = 0; i < fsNeg.size(); ++i) {
                ArrayXXs f = ArrayXXs::Zero(x.size(), y.size());
                cout << fsPos[i].rows() << " =?= " << x.size() << endl;
                cout << fsPos[i].cols() << " =?= " << yPos.size() << endl;
                cout << fsNeg[i].rows() << " =?= " << x.size() << endl;
                cout << fsNeg[i].cols() << " =?= " << yNeg.size() << endl;
                cout << "---" << endl;
                //f.rightCols(yPos.size()) = SQRT1_2 * fsPos[i];
                //f.leftCols(yNeg.size()) = SQRT1_2 * (even ? fsNeg[i] : -fsNeg[i]).rowwise().reverse();
                result.push_back(f);
            }

        }
    }
    return result;
}

template<typename Scalar>
vector<Scalar> SE2DHalf<Scalar>::computeEigenvaluesByIndex(int, int) {
    throw std::logic_error("Function not yet implemented");
}
// std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

template<typename Scalar>
vector<function<Scalar(Scalar, Scalar)>> SE2DHalf<Scalar>::eigenfunctionCalculator(const Scalar &E) {
    vector<function<Scalar(Scalar, Scalar)>> result;
    for (bool even : {false, true}) {
        setParity(even);
        if (abs(E - se2d->findEigenvalue(E)) < 1e-5) {
            for (function<Scalar(Scalar, Scalar)> f : se2d->eigenfunctionCalculator(E)) {
                result.push_back([f, even](const Scalar &x, const Scalar &y) -> Scalar {
                    Scalar r = f(x, y < 0 ? -y : y);
                    if (y < 0 && !even)
                        r = -r;
                    r *= sqrt(Scalar(.5));
                    return r;
                });
            }
        }
    }
    return result;
}


#include "../util/instantiate.h"