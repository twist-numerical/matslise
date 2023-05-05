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
        Scalar error = arg - m_target;
        return -m_margin <= error && error <= m_margin;
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
    Scalar m_target;
    Scalar m_epsilon;
public:
    WithinRel(Scalar target, Scalar epsilon) : m_target(target), m_epsilon(epsilon) {
    }

    bool match(const Scalar &arg) const override {
        if (m_target > m_epsilon || m_target < -m_epsilon) {
            Scalar rel_error = (arg - m_target) / m_target;
            return -m_epsilon <= rel_error && rel_error <= m_epsilon;
        }
        Scalar error = arg - m_target;
        return -m_epsilon <= error && error <= m_epsilon;
    }

    std::string describe() const override {
        Catch::ReusableStringStream sstr;
        sstr << "and " << m_target << " are within " << m_epsilon * 100. << "% of each other";
        return sstr.str();
    }
};


#endif //MATSLISE_TEST_H
