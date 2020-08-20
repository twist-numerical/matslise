#include "../matslise.h"
#include <map>
#include "./matching.h"

using namespace Eigen;
using namespace matslise;
using namespace std;

template<typename Scalar, bool add>
void cec_cce(Matrix<Scalar, Dynamic, Dynamic> &addTo, const Y<Scalar, Dynamic, Dynamic> &y) {
    if constexpr(add)
        addTo += (y).getdY(0).transpose() * (y).getY(1) - (y).getdY(1).transpose() * (y).getY(0);
    else
        addTo -= (y).getdY(0).transpose() * (y).getY(1) - (y).getdY(1).transpose() * (y).getY(0);
}

template<typename Scalar>
vector<Y<Scalar, Dynamic>>
Matslise2D<Scalar>::eigenfunctionSteps(const Y<Scalar, Dynamic> &yLeft, const Scalar &E) const {
    MATSLISE_SCOPED_TIMER("2D eigenfunctionSteps");
    Index sectorCount = sectors.size();
    Index N = config.basisSize;
    auto *steps = new Y<Scalar, Dynamic>[sectorCount + 1];
    auto *endSteps = new Y<Scalar, Dynamic>[sectorCount];

    steps[0] = yLeft;
    steps[sectorCount] = dirichletBoundary;
    auto *U = new MatrixXs[sectorCount];

    for (int i = sectorCount - 1; i > matchIndex; --i) {
        endSteps[i] = i < sectorCount - 1 ? (MatrixXs)(M[i].transpose()) * steps[i + 1] : steps[i + 1];
        if (i + 1 < sectorCount)
            U[i + 1] = conditionY(endSteps[i]);
        steps[i] = sectors[i]->propagate(E, endSteps[i], sectors[i]->max, sectors[i]->min, true);
    }
    const Y<Scalar, Dynamic> matchRight = steps[matchIndex + 1];

    U[0] = MatrixXs::Identity(N, N);
    for (int i = 0; i <= matchIndex; ++i) {
        endSteps[i] = sectors[i]->propagate(E, steps[i], sectors[i]->min, sectors[i]->max, true);
        steps[i + 1] = M[i] * endSteps[i];
        U[i + 1] = conditionY(steps[i + 1]);
    }
    const Y<Scalar, Dynamic> matchLeft = steps[matchIndex + 1];
    MatrixXs kernel = getKernel<Scalar>(matchLeft, matchRight, 1e-5);

    vector<Y<Scalar, Dynamic>> elements;
    if (kernel.cols() > 0) {
        elements.resize(sectorCount + 1);
        MatrixXs left = matchLeft.getY(0).colPivHouseholderQr().solve(kernel);
        MatrixXs right = matchRight.getY(0).colPivHouseholderQr().solve(kernel);

        Y<Scalar, Dynamic> elementMatchRight = matchRight * right;

        MatrixXs normalizer = MatrixXs::Zero(left.cols(), left.cols());
        for (int i = matchIndex + 1; i >= 0; --i) {
            elements[i] = steps[i] * left;
            if (i <= matchIndex) {
                cec_cce<Scalar, true>(normalizer, endSteps[i] * left);
                cec_cce<Scalar, false>(normalizer, elements[i]);
            }
            if (i > 0) {
                U[i].template triangularView<Upper>().
                        template solveInPlace<OnTheLeft>(left);
            }
        }

        for (int i = matchIndex + 1; i < sectorCount; ++i) {
            elements[static_cast<size_t>(i + 1)] = endSteps[i] * right;
            cec_cce<Scalar, true>(normalizer, elements[i + 1]);
            cec_cce<Scalar, false>(normalizer, i > matchIndex + 1 ? steps[i] * right : elementMatchRight);
            if (i + 1 < sectorCount) {
                U[i + 1].template triangularView<Upper>().
                        template solveInPlace<OnTheLeft>(right);
            }
        }

        LLT<Ref<MatrixXs>> llt(normalizer); // In place

        for (Index i = 0; i <= sectorCount; ++i) {
            Y<Scalar, Dynamic> &element = elements[static_cast<size_t>(i)];
            llt.matrixL().transpose().template solveInPlace<OnTheRight>(element.y);
            llt.matrixL().transpose().template solveInPlace<OnTheRight>(element.dy);
        }
    }
    delete[] steps;
    delete[] endSteps;
    delete[] U;
    return elements;
}

template<typename Scalar>
Index findSectorIndex(const Matslise2D<Scalar> *matslise2D, const Scalar &y) {
    Eigen::Index a = 0;
    Eigen::Index b = matslise2D->sectors.size();
    while (a + 1 < b) {
        Eigen::Index c = (a + b) / 2;
        if (y < matslise2D->sectors[c]->min)
            b = c;
        else
            a = c;
    }
    if (a > matslise2D->matchIndex && matslise2D->sectors[a]->min == y)
        a -= 1;
    return a;
}

template<typename Scalar>
template<bool withDerivative>
vector<Eigenfunction2D<Scalar, withDerivative>>
Matslise2D<Scalar>::eigenfunction(const Y<Scalar, Dynamic> &left, const Scalar &E) const {
    shared_ptr<vector<Y<Scalar, Dynamic>>> steps
            = make_shared<vector<Y<Scalar, Dynamic>>>(move(eigenfunctionSteps(left, E)));
    vector<Eigenfunction2D<Scalar, withDerivative>> eigenfunctions;
    if (!steps->empty()) {
        Eigen::Index cols = (*steps)[0].getY(0).cols();
        eigenfunctions.reserve(cols);
        for (Eigen::Index column = 0; column < cols; ++column) {
            eigenfunctions.push_back(
                    {
                            [E, steps, column, this](const Scalar &x, const Scalar &y)
                                    -> typename Eigenfunction2D<Scalar, withDerivative>::ScalarReturn {
                                MATSLISE_SCOPED_TIMER("2D eigenfunction scalar");
                                Index sectorIndex = findSectorIndex(this, y);
                                const Sector *sector = sectors[sectorIndex];

                                Y<Scalar, Dynamic, 1> c =
                                        sector->direction == forward
                                        ? sector->propagate(E, (*steps)[sectorIndex].col(column), sector->min, y)
                                        : sector->propagate(E, (*steps)[sectorIndex + 1].col(column), sector->max, y);

                                if constexpr (withDerivative) {
                                    ArrayXs b, b_x;
                                    tie(b, b_x) = sector->template basis<true>(x);
                                    return {
                                            c.getY(0).dot(b.matrix()),
                                            c.getY(0).dot(b_x.matrix()),
                                            c.getY(1).dot(b.matrix())
                                    };
                                } else {
                                    return c.getY(0).dot(sector->template basis<false>(x).matrix());
                                }
                            },
                            [E, steps, column, this](const ArrayXs &x, const ArrayXs &y)
                                    -> typename Eigenfunction2D<Scalar, withDerivative>::ArrayReturn {
                                MATSLISE_SCOPED_TIMER("2D eigenfunction array");
                                typedef typename std::conditional<withDerivative, std::pair<ArrayXXs, ArrayXXs>, ArrayXXs>::type BasisType;
                                map<Index, BasisType> bases;
                                Index sectorIndex = 0;
                                const Sector *sector = nullptr;
                                const BasisType *basis = nullptr;
                                Index nx = x.size(), ny = y.size();

                                typename Eigenfunction2D<Scalar, withDerivative>::ArrayReturn result;
                                if constexpr(withDerivative) {
                                    result = {ArrayXXs::Zero(nx, ny), ArrayXXs::Zero(nx, ny), ArrayXXs::Zero(nx, ny)};
                                } else {
                                    result = ArrayXXs::Zero(nx, ny);
                                }

                                Y<Scalar, Dynamic, 1> c;
                                for (Index i = 0; i < y.size(); ++i) {
                                    Scalar v = y[i];
                                    if (sector == nullptr || !sector->contains(v)) {
                                        sectorIndex = findSectorIndex(this, v);
                                        sector = sectors[sectorIndex];
                                        if (bases.find(sectorIndex) == bases.end()) {
                                            basis = &(bases[sectorIndex] = sector->template basis<withDerivative>(x));
                                        }
                                    }

                                    if (sector->direction == forward) {
                                        c = sector->propagate(
                                                E, (*steps)[static_cast<size_t>(sectorIndex)].col(column),
                                                sector->min, v, true);
                                    } else {
                                        c = sector->propagate(
                                                E, (*steps)[static_cast<size_t>(sectorIndex + 1)].col(column),
                                                sector->max, v, true);
                                    }


                                    if constexpr (withDerivative) {
                                        MatrixXs phi = get<0>(*basis).matrix() * c.getY(0);
                                        MatrixXs phi_x = get<1>(*basis).matrix() * c.getY(0);
                                        MatrixXs phi_y = get<0>(*basis).matrix() * c.getY(1);

                                        get<0>(result).col(i) = phi;
                                        get<1>(result).col(i) = phi_x;
                                        get<2>(result).col(i) = phi_y;
                                    } else {
                                        result.col(i) = (*basis).matrix() * c.getY(0);
                                    }
                                }

                                return result;
                            }
                    }
            );
        }
    }

    return eigenfunctions;
}

#include "../util/instantiate.h"
