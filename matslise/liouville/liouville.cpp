#include <vector>
#include "../liouville.h"
#include "../util/legendre.h"

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
LiouvilleTransformation<Scalar>::Piece::Piece(
        const LiouvilleTransformation<Scalar> &lt, const Rectangle<Scalar, 1> &r_) :
        r(r_),
        x({Scalar(0), Scalar(0)}) {
    static constexpr const int DEGREE = LiouvilleTransformation<Scalar>::DEGREE;

    Legendre<Scalar, DEGREE + 1> legendreP{lt.m_p, r.min(), r.max()};
    Legendre<Scalar, DEGREE + 1> legendreW{lt.m_w, r.min(), r.max()};

    p = legendreP.asPolynomial();
    w = legendreW.asPolynomial();

    r2x = Legendre<Scalar, DEGREE>{
            [&](Scalar rValue) { return sqrt(w(rValue) / p(rValue)); }, 0, 1
    }.asPolynomial().integral();
    r2x *= (r.max() - r.min());
    x.max() = r2x(1);

    pwx = Legendre<Scalar, DEGREE + 3>{
            [&](Scalar xValue) {
                Scalar rValue = x2r_impl(r2x, xValue);
                return p(rValue) * w(rValue);
            }, x.min(), x.max()
    }.asPolynomial();
}

template<typename Scalar>
void setXMin(typename LiouvilleTransformation<Scalar>::Piece &piece, Scalar xMin) {
    if (xMin != piece.x.min()) {
        Scalar dx = xMin - piece.x.min();
        piece.r2x += dx;
        piece.x.max() += dx;
        piece.x.min() = xMin;
    }
}

template<typename Scalar>
inline Scalar relError(const Scalar &exact, const Scalar &estimate) {
    Scalar err = exact - estimate;
    if (exact > 1 || exact < -1)
        err /= exact;
    return abs(err);
}

template<typename Scalar, typename F, typename Piece=typename LiouvilleTransformation<Scalar>::Piece>
inline Scalar pieceError(const F &f, const Piece &full, const Piece &left, const Piece &right) {
    return relError(f(left) + f(right), f(full));
}

template<typename Scalar>
std::vector<typename LiouvilleTransformation<Scalar>::Piece>
constructPiecewise(const LiouvilleTransformation<Scalar> &lt) {
    using Piece = typename LiouvilleTransformation<Scalar>::Piece;
    std::vector<Piece> pieces;
    std::vector<Piece> toDo;
    toDo.emplace_back(lt, lt.rDomain());
    Scalar xMin = 0;

    while (!toDo.empty()) {
        Piece &next = toDo.back();
        setXMin(next, xMin);

        Scalar rMid = (next.r.min() + next.r.max()) / 2;

        Piece left{lt, {next.r.min(), rMid}};
        setXMin(left, xMin);
        Piece right{lt, {rMid, next.r.max()}};
        setXMin(right, left.x.max());

        if ((next.x.max() - next.x.min()) < 1e-8 || (
                relError(right.x.max(), next.x.max()) < 1e-14
                && pieceError<Scalar>([](const Piece &p) {
                    return p.r2x.integral(1) * p.r.diameter();
                }, next, left, right) < 1e-14
                && pieceError<Scalar>([](const Piece &p) {
                    return Legendre<Scalar, 12>([&p](Scalar x) { return p.x2r(x); }, 0, 1).integrate() * p.x.diameter();
                }, next, left, right) < 1e-14
                && pieceError<Scalar>([](const Piece &p) {
                    return p.pwx.integral(1) * p.x.diameter();
                }, next, left, right) < 1e-14
                && pieceError<Scalar>([](const Piece &p) {
                    return p.pwx(1) - p.pwx(0);
                }, next, left, right) < 1e-12
        )) {
            toDo.pop_back();
            pieces.push_back(left);
            pieces.push_back(right);
            xMin = right.x.max();
        } else {
            toDo.pop_back();
            toDo.push_back(right);
            toDo.push_back(left);
        }
    }

    return pieces;
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::Piece::x2r(Scalar x_) const {
    return r.min() + (r.max() - r.min()) * x2r_impl(r2x, x.min() + x_ * (x.max() - x.min()));
}

template<typename Scalar>
LiouvilleTransformation<Scalar>::LiouvilleTransformation(
        const Rectangle<Scalar, 1> &domain, const Function &p, const Function &q, const Function &w)
        : m_rDomain(domain), m_p(p), m_q(q), m_w(w), pieces(constructPiecewise(*this)) {
}

template<typename Scalar>
const typename LiouvilleTransformation<Scalar>::Piece &
findPiece(const LiouvilleTransformation<Scalar> &lt, Scalar value,
          Rectangle<Scalar, 1> (LiouvilleTransformation<Scalar>::Piece::*f)) {
    auto piece = std::lower_bound(
            lt.pieces.begin(), lt.pieces.end(), value,
            [&f](const typename LiouvilleTransformation<Scalar>::Piece &piece, const Scalar &v) {
                return (piece.*f).max() <= v;
            });
    if (piece == lt.pieces.end()) --piece;

    assert(((*piece).*f).contains(value));

    return *piece;
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::r2x(Scalar r) const {
    const Piece &piece = findPiece(*this, r, &Piece::r);
    return piece.r2x(piece.normalizeR(r));
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::x2r(Scalar x) const {
    const Piece &piece = findPiece(*this, x, &Piece::x);
    return piece.x2r(piece.normalizeX(x));
}

template<typename Scalar>
inline Scalar square(const Scalar &x) {
    return x * x;
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::V(Scalar x) const {
    const Piece &piece = findPiece(*this, x, &Piece::x);

    Scalar nx = piece.normalizeX(x);
    Scalar r = piece.x2r(nx);
    Scalar nr = piece.normalizeR(r);

    Scalar h = piece.x.max() - piece.x.min();
    Scalar pw = piece.pwx(nx);
    Scalar pw_dx = piece.pwx.derivative(nx) / h;
    Scalar pw_ddx = piece.pwx.template derivative<2>(nx) / (h * h);

    return m_q(r) / piece.w(nr) - 3 * square(pw_dx / (4 * pw)) + pw_ddx / (4 * pw);
}

template<typename Scalar>
Y<Scalar> LiouvilleTransformation<Scalar>::z2y(const Scalar &r, const Y<Scalar> &z) const {
    Y<Scalar> y(z);

    const Piece &piece = findPiece(*this, r, &Piece::r);
    Scalar nr = piece.normalizeR(r);
    Scalar p = piece.p(nr);
    Scalar w = piece.w(nr);
    Scalar h = piece.r.diameter();
    Scalar dp = piece.p.derivative(nr) / h;
    Scalar dw = piece.w.derivative(nr) / h;
    Scalar xdr = sqrt(w / p);
    Scalar sigma = 1 / sqrt(p * xdr);
    Scalar dsigma = Scalar(-0.25) * sigma * (dp / p + dw / w);

    y.block(None) *= 1 / sigma;
    y.block(dE) *= 1 / sigma;

    y.block(dX) -= dsigma * y.block(None);
    y.block(dXdE) -= dsigma * y.block(dE);
    y.block(dX) *= 1 / (sigma * xdr);
    y.block(dXdE) *= 1 / (sigma * xdr);

    return y;
}

template<typename Scalar>
Y<Scalar> LiouvilleTransformation<Scalar>::y2z(const Scalar &x, const Y<Scalar> &y) const {
    Y<Scalar> z(y);

    const Piece &piece = findPiece(*this, x, &Piece::x);
    Scalar nr = piece.normalizeR(piece.x2r(piece.normalizeX(x)));
    Scalar p = piece.p(nr);
    Scalar w = piece.w(nr);
    Scalar h = piece.r.diameter();
    Scalar dp = piece.p.derivative(nr) / h;
    Scalar dw = piece.w.derivative(nr) / h;
    Scalar xdr = sqrt(w / p);
    Scalar sigma = 1 / sqrt(p * xdr);
    Scalar dsigma = Scalar(-0.25) * sigma * (dp / p + dw / w);

    z.block(None) *= sigma;
    z.block(dE) *= sigma;

    z.block(dX) *= dsigma;
    z.block(dXdE) *= dsigma;
    z.block(dX) += sigma * xdr * y.block(None);
    z.block(dXdE) += sigma * xdr * y.block(dE);

    return z;
}

#include "./instantiate.h"
