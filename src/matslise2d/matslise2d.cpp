#include <iostream>
#include <map>
#include <queue>
#include "../matslise.h"
#include "../util/quadrature.h"
#include "../util/find_sector.h"
#include "../util/scoped_timer.h"
#include "../util/matching.h"

using namespace Eigen;
using namespace matslise;
using namespace std;
using namespace quadrature;


template<typename Scalar>
Array<Scalar, Dynamic, 1> getGrid(const Scalar &min, const Scalar &max, int count) {
    Array<Scalar, Dynamic, 1> points(count);
    for (int i = 0; i < count; ++i)
        points[i] = min + (max - min) * i / (count - 1);
    return points;
}

template<typename Scalar>
Matslise2D<Scalar>::Matslise2D(const function<Scalar(Scalar, Scalar)> &potential,
                               const matslise::Rectangle<Scalar, 2> &domain, const Config &_config):
        AbstractMatslise2D<Scalar>(potential, domain), MatsliseND<Scalar, Sector>(_config.basisSize), config(_config) {
    MATSLISE_SCOPED_TIMER("2D constructor");
    auto sectorsBuild = sector_builder::getOrAutomatic<Matslise2D<Scalar>, false>(
            _config.ySectorBuilder, _config.tolerance)(this, domain.template min<1>(), domain.template max<1>());
    sectors = std::move(sectorsBuild.sectors);
    matchIndex = sectorsBuild.matchIndex;
    Index sectorCount = sectors.size();
    for (auto &sector : sectors)
        sector->quadratures.reset();

    M.reserve(sectorCount - 1);

    map<pair<int, int>, ArrayXXs> prev;
    map<pair<int, int>, ArrayXXs> next;
    Scalar xWidth = domain.template max<0>() - domain.template min<0>();
    for (int k = 0; k < sectorCount - 1; ++k) {
        if (k > 0) {
            prev = std::move(next);
            next.clear();
        }

        M.push_back(move(
                gauss_kronrod::adaptive<Scalar, MatrixXs, true>([&, k](const ArrayXs &x) {
                    int depth = static_cast<int>(round(log2(xWidth / (x[x.size() - 1] - x[0]))));
                    int offset = static_cast<int>(round(
                            (x[0] - domain.template min<0>()) / (xWidth / (1 << depth))));
                    pair<int, int> key{depth, offset};
                    if (prev.find(key) == prev.end())
                        prev[key] = sectors[k]->template basis<false>(x);
                    const ArrayXXs &prevBasis = prev[key];
                    const ArrayXXs &nextBasis = next[key] = sectors[k + 1]->template basis<false>(x);

                    Array<MatrixXs, Dynamic, 1> result(x.size());
                    for (Index i = 0; i < x.size(); ++i) {
                        result(i) = nextBasis.row(i).matrix().transpose() * prevBasis.row(i).matrix();
                    }
                    return result;
                }, domain.template min<0>(), domain.template max<0>(), 1e-8, [](const MatrixXs &v) {
                    return v.array().abs().maxCoeff();
                })
        ));
    }
}

template<typename Scalar>
Matslise2D<Scalar>::~Matslise2D() {
    for (auto &sector : sectors)
        delete sector;
}


template<typename Scalar>
Scalar Matslise2D<Scalar>::estimatePotentialMinimum() const {
    auto iterator = this->sectors.begin();
    Scalar minimal = (*iterator++)->matslise->estimatePotentialMinimum();
    for (; iterator != this->sectors.end(); ++iterator)
        minimal = min(minimal, (*iterator)->matslise->estimatePotentialMinimum());
    return minimal;
}


#include "../util/instantiate.h"
