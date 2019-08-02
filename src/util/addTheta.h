

#ifndef MATSLISE_ADDTHETA_H
#define MATSLISE_ADDTHETA_H

namespace matslise {
    template<typename Left, typename Scalar>
    Left addTheta(const std::pair <Left, Scalar> &p, Scalar &theta) {
        theta += p.second;
        return p.first;
    }
}

#endif //MATSLISE_ADDTHETA_H
