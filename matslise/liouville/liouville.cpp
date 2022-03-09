#include <iostream>
#include "../liouville.h"
#include "../util/legendre.h"
#include "../util/horner.h"

using namespace matslise;
using namespace matslise::legendre;


template<typename Scalar>
Scalar x2r_impl(const decltype(LiouvilleTransformation<Scalar>::Piece::r2x) &r2x, Scalar x_) {
    Scalar low = 0, high = 1;
    while (high - low > 1e-15) {
        Scalar mid = (high + low) / 2;
        Scalar fMid = r2x(mid);
        if (fMid > x_) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return (low + high) / 2;
}

template<typename Scalar>
struct SimplePiece {
    static constexpr const int DEGREE = LiouvilleTransformation<Scalar>::DEGREE;
    const Scalar min;
    const Scalar max;

    const Legendre<Scalar, DEGREE + 1> p;
    const Legendre<Scalar, DEGREE + 1> w;

    SimplePiece(const LiouvilleTransformation<Scalar> &lt, Scalar rMin, Scalar rMax)
            : min(rMin), max(rMax), p(lt.p, rMin, rMax), w(lt.w, rMin, rMax) {
    }

private:
    mutable std::optional<Polynomial<Scalar, DEGREE >> r2x_cache;
public:
    const Polynomial<Scalar, DEGREE> &r2x() const {
        if (r2x_cache) return *r2x_cache;
        const auto &pp = p.asPolynomial();
        const auto &pw = w.asPolynomial();
        auto &r2x = r2x_cache.emplace(Legendre<Scalar, DEGREE>{
                [&pp, &pw](Scalar r) { return sqrt(pw(r) / pp(r)); }, 0, 1
        }.asPolynomial().integral());
        r2x *= (max - min);
        return r2x;
    }

private:
    mutable std::optional<Polynomial<Scalar, DEGREE >> x2r_cache;
public:
    const Polynomial<Scalar, DEGREE> &x2r() const {
        if (x2r_cache) return *x2r_cache;

        const Polynomial<Scalar, DEGREE> &r2x_p = r2x();

        return x2r_cache.emplace(Legendre<Scalar, DEGREE + 1>{
                [&r2x_p](Scalar x) {
                    return x2r_impl(r2x_p, x);
                }, 0, r2x_p(1)
        }.asPolynomial());
    }
};

template<typename Scalar>
Scalar pwIntegral(const SimplePiece<Scalar> &piece) {
    const constexpr int DEGREE = LiouvilleTransformation<Scalar>::DEGREE;
    const auto &p = piece.p.asPolynomial();
    const auto &w = piece.w.asPolynomial();

    const Polynomial<Scalar, DEGREE> &r2x = piece.r2x();

    return (Legendre<Scalar, DEGREE>{
            [&p, &w, &r2x](Scalar x) {
                Scalar r = x2r_impl(r2x, x);
                Scalar pw = p(r) * w(r);
                return abs(pw);
            }, 0, r2x(1)
    }.integrate());
}

template<typename Scalar>
Scalar r2xIntegral(const SimplePiece<Scalar> &piece, Scalar shift = 0) {
    return (piece.r2x().integral(1) + shift) * (piece.max - piece.min);
}

template<typename Scalar>
Scalar x2rIntegral(const SimplePiece<Scalar> &piece, Scalar shift = 0) {
    return (piece.x2r().integral(1) * (piece.max - piece.min) + shift) * (piece.r2x()(1) - piece.r2x()(0));
}


template<typename Scalar>
Scalar relError(const Scalar &exact, const Scalar &estimate) {
    Scalar err = exact - estimate;
    if (exact > 1 || exact < -1)
        err /= exact;
    return abs(err);
}

template<typename Scalar>
void constructPiecewise(LiouvilleTransformation<Scalar> &lt, const SimplePiece<Scalar> &sp) {
    if ((sp.max - sp.min) > 1e-8) {
        Scalar rMid = (sp.min + sp.max) / 2;

        SimplePiece left{lt, sp.min, rMid};
        SimplePiece right{lt, rMid, sp.max};

        Scalar pBest = left.p.integrate() + right.p.integrate();
        Scalar wBest = left.w.integrate() + right.w.integrate();

        if (
                relError(pBest, sp.p.integrate()) > 1e-15
                || relError(wBest, sp.w.integrate()) > 1e-15
                // || relError(r2xIntegral(left) + r2xIntegral(right, left.r2x()(1)), r2xIntegral(sp)) > 1e-14
                // || relError(x2rIntegral(left) + x2rIntegral(right, left.max - left.min), x2rIntegral(sp)) > 1e-14
                || relError(pwIntegral(left) + pwIntegral(right), pwIntegral(sp)) > 1e-14
                ) {
            // subdivide
            constructPiecewise(lt, left);
            constructPiecewise(lt, right);
            return;
        }
    }

    const constexpr int DEGREE = LiouvilleTransformation<Scalar>::DEGREE;
    Scalar xMin = lt.pieces.empty() ? 0 : lt.pieces.back().x.max;

    Polynomial<Scalar, DEGREE> p = sp.p.asPolynomial();
    Polynomial<Scalar, DEGREE> w = sp.w.asPolynomial();
    Polynomial<Scalar, DEGREE> r2x = sp.r2x();
    Polynomial<Scalar, DEGREE + 2> pwx = Legendre<Scalar, DEGREE + 3>{
            [&p, &w, &r2x](Scalar x) {
                Scalar r = x2r_impl(r2x, x);
                return p(r) * w(r);
            }, 0, r2x(1)
    }.asPolynomial();
    r2x += Polynomial<Scalar, 0>{xMin};
    Scalar xMax = r2x(1);

    lt.pieces.push_back(
            typename LiouvilleTransformation<Scalar>::Piece{
                    .r = {sp.min, sp.max},
                    .x = {xMin, xMax},
                    .r2x = r2x,
                    .p = p,
                    .w = w,
                    .pwx =pwx,
            });
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::Piece::x2r(Scalar x_) const {
    return r.min + (r.max - r.min) * x2r_impl(r2x, x.min + x_ * (x.max - x.min));
}

template<typename Scalar>
LiouvilleTransformation<Scalar>::LiouvilleTransformation(
        const Rectangle<Scalar, 1> &domain, const Function &p, const Function &q, const Function &w)
        : domain(domain), p(p), q(q), w(w) {

    constructPiecewise(*this, SimplePiece{*this, domain.min(), domain.max()});
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::r2x(Scalar r) const {
    auto piece = std::lower_bound(pieces.begin(), pieces.end(), r,
                                  [](const Piece &piece, const Scalar &r) { return piece.r.max <= r; });
    if (piece == pieces.end()) --piece;

    assert(piece->r.min <= r && r <= piece->r.max);

    return piece->r2x(piece->normalizeR(r));
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::x2r(Scalar x) const {
    auto piece = std::lower_bound(pieces.begin(), pieces.end(), x,
                                  [](const Piece &piece, const Scalar &x) { return piece.x.max <= x; });
    if (piece == pieces.end()) --piece;

    assert(piece->x.min <= x && x <= piece->x.max);

    return piece->x2r(piece->normalizeX(x));
}

template<typename Scalar>
inline Scalar square(const Scalar &x) {
    return x * x;
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::V(Scalar x) const {
    auto piece = std::lower_bound(pieces.begin(), pieces.end(), x,
                                  [](const Piece &piece, const Scalar &x) { return piece.x.max <= x; });
    if (piece == pieces.end()) --piece;

    assert(piece->x.min <= x && x <= piece->x.max);

    Scalar nx = piece->normalizeX(x);
    Scalar r = piece->x2r(nx);
    Scalar nr = piece->normalizeR(r);

    Scalar h = piece->x.max - piece->x.min;
    Scalar pw = piece->pwx(nx);
    Scalar pw_dx = piece->pwx.derivative(nx) / h;
    Scalar pw_ddx = piece->pwx.template derivative<2>(nx) / (h * h);

    return q(r) / piece->w(nr) - 3 * square(pw_dx / (4 * pw)) + pw_ddx / (4 * pw);
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::pwx(Scalar x) const {
    auto piece = std::lower_bound(pieces.begin(), pieces.end(), x,
                                  [](const Piece &piece, const Scalar &x) { return piece.x.max <= x; });
    if (piece == pieces.end()) --piece;

    assert(piece->x.min <= x && x <= piece->x.max);
    return piece->pwx(piece->normalizeX(x));
}

#include "./instantiate.h"
