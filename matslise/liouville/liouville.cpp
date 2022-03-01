#include <iostream>
#include "../liouville.h"
#include "../util/legendre.h"
#include "../util/horner.h"

using namespace matslise;
using namespace matslise::legendre;

template<typename Scalar>
struct SimplePiece {
    Scalar min;
    Scalar max;

    Legendre<Scalar, LiouvilleTransformation<Scalar>::DEGREE> p;
    Legendre<Scalar, LiouvilleTransformation<Scalar>::DEGREE> w;

    SimplePiece(const LiouvilleTransformation<Scalar> &lt, Scalar rMin, Scalar rMax)
            : min(rMin), max(rMax), p(lt.p, rMin, rMax), w(lt.w, rMin, rMax) {
    }
};

template<typename Scalar>
Scalar x2r_impl(const decltype(LiouvilleTransformation<Scalar>::Piece::r2x) &r2x, Scalar x_) {
    Scalar low = 0, high = 1;
    while (high - low > 1e-8) {
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
void constructPiecewise(LiouvilleTransformation<Scalar> &lt, const SimplePiece<Scalar> &sp) {
    if ((sp.max - sp.min) > 1e-8) {
        Scalar rMid = (sp.min + sp.max) / 2;

        SimplePiece left{lt, sp.min, rMid};
        SimplePiece right{lt, rMid, sp.max};

        Scalar pBest = left.p.integrate() + right.p.integrate();
        Scalar wBest = left.w.integrate() + right.w.integrate();

        if (std::abs(pBest - sp.p.integrate()) > 1e-13 || std::abs(wBest - sp.w.integrate()) > 1e-13) {
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
    Polynomial<Scalar, DEGREE> r2x = Legendre<Scalar, DEGREE - 1>{
            [&p, &w](Scalar r) { return std::sqrt(w(r) / p(r)); }, 0, 1
    }.asPolynomial().integral();
    r2x *= (sp.max - sp.min);
    r2x += Polynomial<Scalar, 0>{xMin};
    Scalar xMax = r2x(1);

    lt.pieces.push_back(
            typename LiouvilleTransformation<Scalar>::Piece{
                    .r = {sp.min, sp.max},
                    .x = {xMin, xMax},
                    .r2x = r2x,
                    .p = p,
                    .w = w,
                    .pwx = Legendre<Scalar, DEGREE>{
                            [&p, &w, &r2x](Scalar x) {
                                Scalar r = x2r_impl(r2x, x);
                                return p(r) * w(r);
                            }, xMin, xMax
                    }.asPolynomial(),
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

    return piece->r2x(piece->normalizeR(r));
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::x2r(Scalar x) const {
    auto piece = std::lower_bound(pieces.begin(), pieces.end(), x,
                                  [](const Piece &piece, const Scalar &x) { return piece.x.max <= x; });
    if (piece == pieces.end()) --piece;

    return piece->x2r(piece->normalizeX(x));
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::V(Scalar x) const {
    auto piece = std::lower_bound(pieces.begin(), pieces.end(), x,
                                  [](const Piece &piece, const Scalar &x) { return piece.x.max <= x; });
    if (piece == pieces.end()) --piece;

    Scalar nx = piece->normalizeX(x);
    Scalar r = piece->x2r(nx);
    Scalar nr = piece->normalizeR(r);

    Scalar h = piece->x.max - piece->x.min;
    Scalar pw = piece->pwx(nx);
    Scalar pw_dx = piece->pwx.derivative(nx) / h;
    Scalar pw_ddx = piece->pwx.template derivative<2>(nx) / (h * h);

    return q(r) / piece->w(nr) - (3 * pw_dx * pw_dx - 4 * pw * pw_ddx) / (16 * pw * pw);
}

#include "./instantiate.h"
