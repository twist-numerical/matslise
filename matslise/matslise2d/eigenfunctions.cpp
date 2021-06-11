#include "../matslise.h"
#include <map>
#include "../util/matching.h"

using namespace Eigen;
using namespace matslise;
using namespace std;


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
    shared_ptr<vector<Y<Scalar, Dynamic>>> steps = make_shared<vector<Y<Scalar, Dynamic>>>(
            move(MatsliseND<Scalar, Matslise2DSector<Scalar>>::eigenfunctionSteps(left, E)));
    vector<Eigenfunction2D<Scalar, withDerivative>> eigenfunctions;
    if (!steps->empty()) {
        Eigen::Index cols = (*steps)[0].block().cols();
        eigenfunctions.reserve(cols);
        for (Eigen::Index column = 0; column < cols; ++column) {
            eigenfunctions.push_back(
                    {
                            [E, steps, column, this](const Scalar &x, const Scalar &y)
                                    -> typename Eigenfunction2D<Scalar, withDerivative>::ScalarReturn {
                                MATSLISE_SCOPED_TIMER("2D eigenfunction scalar");
                                Index sectorIndex = findSectorIndex(this, y);
                                const value_ptr<Sector> &sector = sectors[sectorIndex];

                                Y<Scalar, Dynamic, 1> c =
                                        sector->direction == forward
                                        ? sector->propagate(E, (*steps)[sectorIndex].col(column), sector->min, y)
                                        : sector->propagate(E, (*steps)[sectorIndex + 1].col(column), sector->max, y);

                                if constexpr (withDerivative) {
                                    ArrayXs b, b_x;
                                    tie(b, b_x) = sector->template basis<true>(x);
                                    return {
                                            c.block().dot(b.matrix()),
                                            c.block().dot(b_x.matrix()),
                                            c.block(YDiff::dX).dot(b.matrix())
                                    };
                                } else {
                                    return c.block().dot(sector->template basis<false>(x).matrix());
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
                                        sector = sectors[sectorIndex].get();
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
                                        MatrixXs phi = get<0>(*basis).matrix() * c.block();
                                        MatrixXs phi_x = get<1>(*basis).matrix() * c.block();
                                        MatrixXs phi_y = get<0>(*basis).matrix() * c.block(YDiff::dX);

                                        get<0>(result).col(i) = phi;
                                        get<1>(result).col(i) = phi_x;
                                        get<2>(result).col(i) = phi_y;
                                    } else {
                                        result.col(i) = (*basis).matrix() * c.block();
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
