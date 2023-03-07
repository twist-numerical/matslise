#include "../test.h"
#include "test_problem.h"
#include "../../matslise/matslise.h"

using matslise::Y;
using matslise::Matslise;
using matslise::MatsliseHalf;

template<typename Scalar>
void testHarmonic(Scalar tolerance) {
    Matslise<Scalar> ms([](Scalar x) -> Scalar {
        return x * x;
    }, -15, 15, 0.001 * tolerance);

    std::vector<Scalar> correct;
    for (int i = 0; i < 50; ++i)
        correct.push_back(2 * i + 1);

    testProblem<Scalar>(ms, Y<Scalar>::Dirichlet(), Y<Scalar>::Dirichlet(), correct, tolerance, 1e-6);
}


TEST_CASE("harmonic oscillator", "[matslise][harmonic]") {
    testHarmonic<double>(1e-7);
}


#ifdef MATSLISE_LONG_DOUBLE

TEST_CASE("harmonic oscillator (long)", "[matslise][harmonic][long]") {
    testHarmonic<long double>(1e-11);
}

#endif

#ifdef MATSLISE_QUADMATH

#include <boost/multiprecision/float128.hpp>

using boost::multiprecision::float128;

TEST_CASE("harmonic oscillator (float128)", "[matslise][harmonic][float128]") {
    testHarmonic<float128>(1e-14);
}

#endif
