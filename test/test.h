#ifndef MATSLISE_TEST_H
#define MATSLISE_TEST_H

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinRel;

template<typename Scalar>
struct WithinAbs : public Catch::Matchers::MatcherBase<Scalar> {
    Scalar m_target;
    Scalar m_margin;

    WithinAbs(Scalar target, Scalar margin) : m_target(target), m_margin(margin) {
    }

    bool match(const Scalar &arg) const override {
        return abs(m_target - arg) <= m_margin;
    }

    std::string describe() const override {
        Catch::ReusableStringStream sstr;
        sstr << "is within " << m_margin << " of " << m_target;
        return sstr.str();
    }
};


#endif //MATSLISE_TEST_H
