#ifndef MATSLISE_LIOUVILLE_H
#define MATSLISE_LIOUVILLE_H

#include <vector>
#include "./util/eigen.h"
#include "./util/rectangle.h"
#include "./util/polynomial.h"

namespace matslise {
    template<typename Scalar>
    class LiouvilleTransformation {
    public:
        static constexpr const int DEGREE = 8;

        typedef std::function<Scalar(Scalar)> Function;

        const Rectangle<Scalar, 1> rDomain;
        const Function p, q, w;

        struct Piece {
            Rectangle<Scalar, 1> r;
            Rectangle<Scalar, 1> x;

            Piece(const LiouvilleTransformation<Scalar> &, const Rectangle<Scalar, 1>&);

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

        Rectangle<Scalar, 1> xDomain() const {
            return {pieces.front().x.min(), pieces.back().x.max()};
        }
    };

}

#endif //MATSLISE_LIOUVILLE_H
