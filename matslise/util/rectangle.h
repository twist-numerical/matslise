#ifndef MATSLISE_RECTANGLE_H
#define MATSLISE_RECTANGLE_H

#include <cassert>
#include <array>

namespace matslise {

    template<typename Scalar, int n>
    struct Rectangle {
    private:
        std::array<Scalar, n * 2> bounds;
    public:
        constexpr Rectangle() {};

        template<typename... T, typename=typename std::enable_if<
                sizeof...(T) == 2 * n && (std::is_convertible<T, Scalar>::value && ...)>::type>
        constexpr Rectangle(T... bounds) : bounds({bounds...}) {}

        constexpr Scalar min(int axis) const {
            return bounds[2 * axis];
        }

        Scalar &min(int axis) {
            return bounds[2 * axis];
        }

        template<int axis = n == 1 ? 0 : -1>
        constexpr Scalar min() const {
            static_assert(0 <= axis && axis < n, "Specify a valid axis");
            return bounds[2 * axis];
        }

        template<int axis = n == 1 ? 0 : -1>
        Scalar &min() {
            static_assert(0 <= axis && axis < n, "Specify a valid axis");
            return bounds[2 * axis];
        }

        constexpr Scalar max(int axis) const {
            return bounds[2 * axis + 1];
        }

        Scalar &max(int axis) {
            return bounds[2 * axis + 1];
        }

        template<int axis = n == 1 ? 0 : -1>
        constexpr Scalar max() const {
            static_assert(0 <= axis && axis < n, "Specify a valid axis");
            return bounds[2 * axis + 1];
        }

        template<int axis = n == 1 ? 0 : -1>
        Scalar &max() {
            static_assert(0 <= axis && axis < n, "Specify a valid axis");
            return bounds[2 * axis + 1];
        }

        constexpr bool contains(int axis, const Scalar &v) const {
            return min(axis) <= v && v <= max(axis);
        }

        template<int axis = n == 1 ? 0 : -1>
        constexpr bool contains(const Scalar &v) const {
            static_assert(0 <= axis && axis < n, "Specify a valid axis");
            return min<axis>() <= v && v <= max<axis>();
        }

        constexpr bool contains(const std::array<Scalar, n> &point) const {
            for (int i = 0; i < n; ++i)
                if (!contains(i, point[i]))
                    return false;
            return true;
        }

        Scalar diameter() const {
            if constexpr(n == 1)
                return abs(max<0>() - min<0>());
            if constexpr(n == 2)
                return hypot(max<0>() - min<0>(), max<1>() - min<1>());
            Scalar s = 0;
            for (int i = 0; i < n; ++i) {
                Scalar d = max(i) - min(i);
                s += d * d;
            }
            return sqrt(s);
        }

        template<typename... T>
        constexpr typename std::enable_if<(std::is_integral<T>::value && ...), Rectangle<Scalar, sizeof...(T)>>::type
        slice(T... axes) const {
            Rectangle<Scalar, sizeof...(T)> r;
            int j = 0;
            (..., (r.min(j) = min(axes), r.max(j) = max(axes), ++j));
            return r;
        }

        template<int... axes>
        constexpr Rectangle<Scalar, sizeof...(axes)> slice() const {
            Rectangle<Scalar, sizeof...(axes)> r;
            int j = 0;
            (..., (r.min(j) = min(axes), r.max(j) = max(axes), ++j));
            return r;
        }

        constexpr bool operator==(const Rectangle<Scalar, n> &rhs) const {
            return bounds == rhs.bounds;
        }

        constexpr bool operator!=(const Rectangle<Scalar, n> &rhs) const {
            return !(*this == rhs);
        }
    };
}

#endif //MATSLISE_RECTANGLE_H
