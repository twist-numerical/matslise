#ifndef MATSLISE_TEST_H
#define MATSLISE_TEST_H

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

template<typename Scalar>
class WithinAbs : public Catch::Matchers::MatcherBase<Scalar> {
private:
    Scalar m_target;
    Scalar m_margin;

public:
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

template<typename Scalar>
class WithinRel : public Catch::Matchers::MatcherBase<Scalar> {
private:
    double m_target;
    double m_epsilon;
public:
    WithinRel(Scalar target, Scalar epsilon) : m_target(target), m_epsilon(epsilon) {
    }

    bool match(const Scalar &arg) const override {
        if (m_target > m_epsilon || m_target < -m_epsilon)
            return abs((arg - m_target) / m_target) <= m_epsilon;
        return abs(arg) <= m_epsilon;
    }

    std::string describe() const override {
        Catch::ReusableStringStream sstr;
        sstr << "and " << m_target << " are within " << m_epsilon * 100. << "% of each other";
        return sstr.str();
    }
};


#endif //MATSLISE_TEST_H
