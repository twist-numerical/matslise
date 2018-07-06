//
// Created by toon on 6/28/18.
//

#ifndef MATSLISE_MATRIX2D_H
#define MATSLISE_MATRIX2D_H

#include <ostream>
#include <Eigen/Dense>

template<typename V = double>
class Vector2D {
public:
    V x, y;

    Vector2D() {
    };

    Vector2D(std::initializer_list<V> xy) {
        auto i = xy.begin();
        x = *i;
        y = *(++i);
    };

    Vector2D(const Eigen::Matrix<V, 2, 1> &r) {
        x = r(0);
        y = r(1);
    };

    Vector2D<V> operator-() const {
        return {-x, -y};
    }

    Vector2D<V> operator+(const Vector2D<V> r) const {
        return {x + r.x, y + r.y};
    }

    Vector2D<V> operator-(const Vector2D<V> r) const {
        return {x - r.x, y - r.y};
    }

    Vector2D<V> operator*(const V &f) const {
        return {x * f, y * f};
    }

    Vector2D<V> operator*=(const V &f) {
        x *= f;
        y *= f;
        return *this;
    }

    Vector2D<V> operator/(const V &f) const {
        return {x / f, y / f};
    }

    Vector2D<V> &operator=(const Eigen::Matrix<V, 2, 1> &r) {
        x = r(0);
        y = r(1);
        return *this;
    }

    V operator[](int a) const {
        return a == 0 ? x : y;
    }

    operator Eigen::Matrix<V, 2, 1>() {
        Eigen::Matrix<V, 2, 1> e;
        e << x, y;
        return e;
    };

    friend Vector2D<V> operator*(const V &f, const Vector2D<V> &v) {
        return {f * v.x, f * v.y};
    }
};

template<typename V = double>
class Matrix2D {
public:
    V a, b, c, d;

    Matrix2D<V> operator+(const Matrix2D<V> &r) const {
        return {a + r.a, b + r.b, c + r.c, d + r.d};
    }

    Matrix2D<V> operator+=(const Matrix2D<V> &r) {
        a += r.a;
        b += r.b;
        c += r.c;
        d += r.d;
        return *this;
    }

    Matrix2D<V> operator-(const Matrix2D<V> &r) const {
        return {a - r.a, b - r.b, c - r.c, d - r.d};
    }

    Matrix2D<V> operator*(const Matrix2D<V> &r) const {
        return {a * r.a + b * r.c, a * r.b + b * r.d, c * r.a + d * r.c, c * r.b + d * r.d};
    }

    Matrix2D<V> operator*(const V &f) const {
        return {a * f, b * f, c * f, d * f};
    }

    Vector2D<V> operator*(const Vector2D<V> &r) const {
        return {a * r.x + b * r.y, c * r.x + d * r.y};
    }

    Matrix2D<V> operator*=(double f) {
        a *= f;
        b *= f;
        c *= f;
        d *= f;
        return *this;
    }


    friend Matrix2D<V> operator*(const V &f, const Matrix2D<V> &m) {
        return {f * m.a, f * m.b, f * m.c, f * m.d};
    }

    friend std::ostream &operator<<(std::ostream &os, const Matrix2D<V> &m) {
        return os << "/" << m.a << ", " << m.b << "\\\n\\" << m.c << ", " << m.d << "/";
    }
};

#endif //MATSLISE_MATRIX2D_H
