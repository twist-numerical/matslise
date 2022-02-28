#ifndef MATSLISE_LIOUVILLE_H
#define MATSLISE_LIOUVILLE_H

#include "./util/eigen.h"
#include "./util/rectangle.h"

namespace matslise {
    template<typename Scalar>
    class LiouvilleTransformation {
    public:
        static constexpr const int legendrePoints = 6;

        typedef std::function<Scalar(Scalar)> Function;

        const Rectangle<Scalar, 1> domain;
        const Function p, q, w;

        struct RXPoint {
            struct {
                Scalar min;
                Scalar max;
            } r;
            struct {
                Scalar min;
                Scalar max;
            } x;
            std::array<Scalar, legendrePoints> legendre;

            Scalar r2x(Scalar r_) const;

            Scalar x2r(Scalar x_) const;
        };

        const std::vector<RXPoint> rxMap;

        LiouvilleTransformation(const Rectangle<Scalar, 1> &domain, const Function &p, const Function &q,
                                const Function &w);

        Scalar r2x(Scalar r) const;

        Scalar x2r(Scalar x) const;
    };

}

#endif //MATSLISE_LIOUVILLE_H
