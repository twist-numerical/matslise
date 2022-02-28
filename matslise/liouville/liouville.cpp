#include "../liouville.h"
#include "../util/legendre.h"
#include "../util/horner.h"

using namespace matslise;
using namespace matslise::legendre;

template<typename Scalar>
struct LT {
    typedef Legendre<LiouvilleTransformation<Scalar>::legendrePoints, Scalar, Scalar> Quadrature;
};

template<typename Scalar, typename F>
typename LT<Scalar>::Quadrature
potentialIntegrate(const LiouvilleTransformation<Scalar> &lt, Scalar min, Scalar max) {
    return typename LT<Scalar>::Quadrature{
            [&lt](Scalar rPrime) {
                return std::sqrt(lt.w(rPrime) * lt.p(rPrime));
            }, min, max};
}

template<typename Scalar>
Scalar constructRXMapping(LiouvilleTransformation<Scalar> &lt, Scalar rMin, Scalar rMax, Scalar xOffset,
                          const typename LT<Scalar>::Quadrature &legendre, const Scalar &integral) {
    Scalar rMid = (rMin + rMax) / 2;
    auto legLeft = potentialIntegrate(lt, rMin, rMid);
    Scalar qLeft = legLeft.integrate();
    auto legRight = potentialIntegrate(lt, rMid, rMax);
    Scalar qRight = legRight.integrate();

    Scalar qBest = qLeft + qRight;
    if (std::abs(qBest - integral) < 1e-10) {
        lt.rxMap.push_back({.r = rMax, .x = xOffset + qBest, .legendre=legendre.asPolynomial()});
        return qBest;
    } else {
        qLeft = constructRXMapping(lt, rMin, rMid, xOffset, legLeft, qLeft);
        qRight = constructRXMapping(lt, rMid, rMax, xOffset + qLeft, legRight, qRight);
        return qLeft + qRight;
    }
}


template<typename Scalar>
LiouvilleTransformation<Scalar>::LiouvilleTransformation(
        const Rectangle<Scalar, 1> &domain, const Function &p, const Function &q, const Function &w)
        : domain(domain), p(p), q(q), w(w) {
    rxMap.push_back({.r=domain.min(), .x = 0});

    auto leg = potentialIntegrate<Scalar>(*this, domain.min(), domain.max());
    constructRXMapping(*this, domain.min(), domain.max(), 0, leg, leg.integrate());
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::r2x(Scalar r) const {
    auto lower = std::lower_bound(rxMap.begin(), rxMap.end(), r,
                                  [](const RXPoint &rx, const Scalar &r) { return rx.r <= r; });
    if (lower == rxMap.end()) --lower;

    return lower->r2x(r);
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::x2r(Scalar x) const {
    auto lower = std::lower_bound(rxMap.begin(), rxMap.end(), x,
                                  [](const RXPoint &rx, const Scalar &x) { return rx.x <= x; });
    if (lower == rxMap.end()) --lower;

    return lower->x2r(x);
}


template<typename Scalar>
Scalar r2x_normalized(const typename LiouvilleTransformation<Scalar>::RXPoint &rx, Scalar nr) {
    const int lp = LiouvilleTransformation<Scalar>::legendrePoints;
    Scalar x_ = rx.legendre[lp - 1] / lp;
    for (int i = lp - 2; i >= 0; --i)
        x_ = nr * x_ + rx.legendre[i] / (i + 1);
    return x_ * nr;
}


template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::RXPoint::r2x(Scalar r_) const {
    Scalar nr = (r_ - r.min) / (r.max - r.min);
    return r2x_normalized(*this, nr);
}

template<typename Scalar>
Scalar LiouvilleTransformation<Scalar>::RXPoint::x2r(Scalar x_) const {
    Scalar low = 0, high = 1;
    while (high - low > 1e-8) {
        Scalar mid = (high + low) / 2;
        Scalar fMid = r2x_normalized(mid);
        if (fMid < x_) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return (low + high) / 2;
}
