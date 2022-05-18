#ifndef MATSLISE_LIOUVILLE_H
#define MATSLISE_LIOUVILLE_H

#include <vector>
#include "./util/eigen.h"
#include "./util/rectangle.h"
#include "./util/polynomial.h"
#include "./matslise.h"

namespace matslise {
    template<typename Scalar>
    class LiouvilleTransformation {
    public:
        using Function = std::function<Scalar(Scalar)>;
    private:
        static constexpr const int DEGREE = 8;

        const Rectangle<Scalar, 1> m_rDomain;
        const Function m_p, m_q, m_w;
    public:

        struct Piece {
            Rectangle<Scalar, 1> r;
            Rectangle<Scalar, 1> x;

            Piece(const LiouvilleTransformation<Scalar> &, const Rectangle<Scalar, 1> &);

            Polynomial<Scalar, DEGREE> r2x; // [0,1] -> [xmin, xmax]
            Polynomial<Scalar, DEGREE> p; // [0,1] -> p(r)
            Polynomial<Scalar, DEGREE> w; // [0,1] -> w(r)
            Polynomial<Scalar, DEGREE + 2> pwx; // [0,1] -> p(x)*w(x)

            Scalar normalizeR(Scalar r_) const {
                return (r_ - r.min()) / (r.max() - r.min());
            };

            Scalar normalizeX(Scalar x_) const {
                return (x_ - x.min()) / (x.max() - x.min());
            };

            Scalar x2r(Scalar x_) const;
        };

        std::vector<Piece> pieces;

        LiouvilleTransformation(const Rectangle<Scalar, 1> &domain, const Function &p, const Function &q,
                                const Function &w);

        Scalar r2x(Scalar r) const;

        Scalar x2r(Scalar x) const;

        Scalar V(Scalar x) const;

        // This does not preserve scaling
        Y<Scalar> transformY(const Scalar &r, const Y<Scalar> &yr) const;

        Rectangle<Scalar, 1> xDomain() const {
            return {pieces.front().x.min(), pieces.back().x.max()};
        }

        const Rectangle<Scalar, 1> &rDomain() const {
            return m_rDomain;
        }
    };

    template<typename Scalar>
    class SturmLiouville {
        using Function = typename LiouvilleTransformation<Scalar>::Function;
        LiouvilleTransformation<Scalar> transformation;
        Matslise<Scalar> matslise;
    public:
        SturmLiouville(const Function &p, const Function &q, const Function &w,
                       const Rectangle<Scalar, 1> &domain, const Scalar &tolerance)
                : transformation(domain, p, q, w),
                  matslise(std::bind(&LiouvilleTransformation<Scalar>::V, transformation, std::placeholders::_1),
                           transformation.xDomain(), tolerance) {
        }

        std::vector<std::pair<int, Scalar>>
        eigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const {
            return matslise.eigenvaluesByIndex(
                    Imin, Imax,
                    transformation.transformY(transformation.rDomain().min(), left),
                    transformation.transformY(transformation.rDomain().max(), right));
        }
    };

}

#endif //MATSLISE_LIOUVILLE_H
