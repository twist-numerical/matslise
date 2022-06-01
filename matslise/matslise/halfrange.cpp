#include <iostream>
#include "../matslise.h"
#include "../util/constants.h"
#include <memory>

using namespace std;
using namespace matslise;
using namespace Eigen;

#define EPS (1e-12)

template<typename Scalar>
MatsliseHalf<Scalar>::MatsliseHalf(
        function<Scalar(Scalar)> V, const Scalar &xmax, const Scalar &tolerance,
        SectorBuilder<Matslise<Scalar>> sectorBuilder
): AbstractMatslise<Scalar>(V, {-xmax, xmax}) {
    ms.reset(new Matslise<Scalar>(V, {Scalar(0), xmax}, tolerance, sectorBuilder));
}

template<typename Scalar>
inline bool isEven(const MatsliseHalf<Scalar> *hr, Scalar E, const Y<Scalar> &side, int index) {
    if (index == -1) {
        Scalar error0 = get<0>(hr->ms->matchingError(E, Y<Scalar>::Neumann(), side));
        Scalar error1 = get<0>(hr->ms->matchingError(E, Y<Scalar>::Dirichlet(), side));
        return abs(error0) < abs(error1);
    }
    return bool(index % 2 == 0);
}

template<typename Scalar>
inline Y<Scalar> getY0(bool even) {
    return even ? Y<Scalar>::Neumann() : Y<Scalar>::Dirichlet();
}


template<typename Scalar>
class EigenfunctionHalf : public AbstractMatslise<Scalar>::Eigenfunction {
public:
    bool even;
    unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction> eigenfunction;

    EigenfunctionHalf(bool even, unique_ptr<typename AbstractMatslise<Scalar>::Eigenfunction> f)
            : even(even), eigenfunction(std::move(f)) {
    }

    virtual Array<Scalar, 2, 1> operator()(const Scalar &x) const override {
        Array<Scalar, 2, 1> c = (*eigenfunction)(x < 0 ? -x : x);
        c *= sqrt(Scalar(.5));
        if (x < 0 && !even)
            c *= -1;
        return c;
    }

    virtual Array<Scalar, Dynamic, 2>
    operator()(const Array<Scalar, Dynamic, 1> &x) const override {
        Eigen::Index n = x.size();
        for (Eigen::Index i = 1; i < n; ++i)
            if (x[i - 1] > x[i])
                throw runtime_error("Matslise::computeEigenfunction(): x has to be sorted");

        Eigen::Index negatives = 0;
        for (Eigen::Index i = 0; i < n; ++i)
            if (x[i] < 0)
                negatives = i + 1;

        Array<Scalar, Dynamic, 1> xNeg(negatives);
        Array<Scalar, Dynamic, 1> xPos(n - negatives);

        for (Eigen::Index i = 0; i < negatives; ++i)
            xNeg[i] = -x[negatives - 1 - i];

        for (Eigen::Index i = negatives; i < n; ++i)
            xPos[i - negatives] = x[i];

        Array<Scalar, Dynamic, 2> yNeg, yPos, ys(n, 2);
        yNeg = (*eigenfunction)(xNeg);
        yPos = (*eigenfunction)(xPos);

        static const Scalar SQRT1_2 = sqrt(Scalar(.5));
        for (Eigen::Index i = 0; i < negatives; ++i) {
            Scalar f = (even ? 1 : -1) * SQRT1_2;
            (ys.row(negatives - 1 - i) = yNeg.row(i) * f)[1] *= -1;
        }
        for (Eigen::Index i = negatives; i < n; ++i)
            ys.row(i) = yPos.row(i - negatives) * SQRT1_2;

        return ys;
    }
};

template<typename Scalar>
typename std::unique_ptr<typename MatsliseHalf<Scalar>::Eigenfunction>
MatsliseHalf<Scalar>::eigenfunction(const Scalar &E, const Y<Scalar> &side, int index) const {
    bool even = isEven(this, E, side, index);
    return std::make_unique<EigenfunctionHalf<Scalar>>(even, ms->eigenfunction(E, getY0<Scalar>(even), side));
}

template<typename Scalar>
vector<pair<int, Scalar>>
mergeEigenvalues(const vector<pair<int, Scalar>> &even, const vector<pair<int, Scalar>> &odd) {
    vector<pair<int, Scalar>> values;

    auto a = even.begin();
    auto b = odd.begin();
    while (a != even.end() || b != odd.end()) {
        if (a == even.end() || (b != odd.end() && b->first < a->first)) {
            values.push_back({2 * b->first + 1, b->second});
            ++b;
        } else {
            values.push_back({2 * a->first, a->second});
            ++a;
        }
    }

    return values;
}

template<typename Scalar>
Scalar MatsliseHalf<Scalar>::eigenvalueError(const Scalar &E, const Y<Scalar> &side, int index) const {
    return ms->eigenvalueError(E, getY0<Scalar>(isEven(this, E, side, index)), side);
}

template<typename Scalar>
vector<pair<int, Scalar>>
MatsliseHalf<Scalar>::eigenvaluesByIndex(int Imin, int Imax, const Y<Scalar> &side) const {
    return mergeEigenvalues(
            ms->eigenvaluesByIndex(Imin / 2 + Imin % 2, Imax / 2 + Imax % 2, Y<Scalar>::Neumann(), side),
            ms->eigenvaluesByIndex(Imin / 2, Imax / 2, Y<Scalar>::Dirichlet(), side));
}

template<typename Scalar>
vector<pair<int, Scalar>>
MatsliseHalf<Scalar>::eigenvalues(
        const Scalar &Emin, const Scalar &Emax, const Y<Scalar> &side) const {
    return mergeEigenvalues(ms->eigenvalues(Emin, Emax, Y<Scalar>::Neumann(), side),
                            ms->eigenvalues(Emin, Emax, Y<Scalar>::Dirichlet(), side));
}

#include "instantiate.h"