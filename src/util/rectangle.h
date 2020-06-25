#ifndef MATSLISE_RECTANGLE_H
#define MATSLISE_RECTANGLE_H

#include <cassert>

namespace matslise {

    template<int n, typename Scalar>
    class Rectangle {
    public:
        Rectangle<n - 1, Scalar> sub;
        Scalar min, max;

        constexpr Scalar getMin(int axis) const {
            if (axis + 1 == n)
                return min;
            return sub.getMin(axis);
        }

        template<int axis>
        constexpr Scalar getMin() const {
            if constexpr (axis + 1 == n)
                return min;
            return sub.template getMin<axis>();
        }

        constexpr Scalar getMax(int axis) const {
            if (axis + 1 == n)
                return max;
            return sub.getMax(axis);
        }

        template<int axis>
        constexpr Scalar getMax() const {
            if constexpr (axis + 1 == n)
                return max;
            return sub.template getMax<axis>();
        }

        constexpr bool contains(int axis, const Scalar &v) const {
            if (axis + 1 == n)
                return min <= v && v <= max;
            return sub.contains(axis, v);
        }

        Scalar diameter() const {
            return hypot(max - min, sub.diameter());
        }
    };

    template<typename Scalar>
    class Rectangle<1, Scalar> {
    public:
        Scalar min, max;

        constexpr Scalar getMin(int axis) const {
            assert(axis == 0);
            (void) (axis); // UNUSED
            return min;
        }

        template<int axis>
        constexpr Scalar getMin() const {
            assert(axis == 0);
            return min;
        }

        constexpr Scalar getMax(int axis) const {
            assert(axis == 0);
            (void) (axis); // UNUSED
            return max;
        }

        template<int axis>
        constexpr Scalar getMax() const {
            assert(axis == 0);
            return max;
        }

        constexpr bool contains(int axis, const Scalar &v) const {
            assert(axis == 0);
            (void) (axis); // UNUSED
            return contains(v);
        }

        constexpr bool contains(const Scalar &v) const {
            return min <= v && v <= max;
        }

        Scalar diameter() const {
            return abs(max - min);
        }
    };
}

#endif //MATSLISE_RECTANGLE_H
