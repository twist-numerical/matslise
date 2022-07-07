#include "../catch.hpp"
#include "../../matslise/matslise.h"
#include <map>
#include <set>

using namespace matslise;
using namespace std;
using namespace Eigen;

/* TEST_CASE("Periodic zero potential", "[matslise][periodic]") {
    PeriodicMatslise<> matslise([](double) { return 0; }, 0, M_PI, 1e-8);
    auto error = matslise.matchingError(3);

}

TEST_CASE("Periodic mathieu", "[matslise][periodic]") {
    PeriodicMatslise<> matslise([](double x) { return 2 * cos(2 * x); }, 0, M_PI, 1e-8);

    for (auto &eigenpair : matslise.eigenpairsByIndex(0, 20)) {
        cout << "i: " << get<0>(eigenpair) << ", E: " << get<1>(eigenpair) << ", m: " << get<2>(eigenpair).size()
             << endl;
    }
} */

TEST_CASE("Periodic Andrews asymmetric", "[matslise][periodic]") {
    // https://www.sciencedirect.com/science/article/abs/pii/0168927495000675
    // https://www.cambridge.org/core/journals/anziam-journal/article/correction-of-finite-difference-eigenvalues-of-periodic-sturmliouville-problems/942FF24F2A3B577649B12D7CD17EFDDE
    std::vector<double> exact{
            2.02942, 6.50049, 7.01506, 18.58477, 18.66548, 38.58162, 38.62154, 66.58204, 66.60537,
            102.58252, 102.59772, 146.58286, 146.59352, 198.58310, 198.59998, 258.58326, 258.58931, 326.58338,
            326.58817, 402.58347
    };

    PeriodicMatslise<> matslise([](double x) { return x * x * (M_PI - x); }, 0, M_PI, 1e-6);

    {
        auto eigenvalues = matslise.eigenvaluesByIndex(0, 20);
        REQUIRE(eigenvalues.size() >= 20); // All eigenvalues should be single, more may be found
        for (int i = 0; i < 20; ++i) {
            auto &eigenvalue = eigenvalues[i];
            REQUIRE(get<0>(eigenvalue) == i);
            REQUIRE(Approx(exact[i]).epsilon(1e-4) == get<1>(eigenvalue));
            REQUIRE(get<2>(eigenvalue) == 1);
        }
    }

    {
        auto eigenvalues = matslise.eigenvalues(10, 102.59);
        REQUIRE(eigenvalues.size() == 7); // All eigenvalues should be single, only the required ones should be found
        auto ei = eigenvalues.begin();
        for (int i = 3; i < 10; ++i, ++ei) {
            REQUIRE(get<0>(*ei) == i);
            REQUIRE(Approx(exact[i]).epsilon(1e-4) == get<1>(*ei));
            REQUIRE(get<2>(*ei) == 1);
        }
    }
}

TEST_CASE("Periodic Andrews symmetric", "[matslise][periodic]") {
    // https://www.cambridge.org/core/journals/anziam-journal/article/correction-of-finite-difference-eigenvalues-of-periodic-sturmliouville-problems/942FF24F2A3B577649B12D7CD17EFDDE
    std::map<int, double> exact{
            {1,  2.09946},
            {2,  7.44911},
            {3,  16.64822},
            {4,  17.09658},
            {5,  36.35887},
            {6,  36.36090},
            {7,  64.19884},
            {9,  100.12637},
            {11, 144.08745},
            {13, 196.06412},
            {15, 256.04903},
            {17, 324.03870},
            {19, 400.03133},
            {29, 900.01390},
            {37, 1444.00866},
            {39, 1600.00782}
    };
    std::set<int> seen;

    PeriodicMatslise<> matslise([](double x) { return 10 * cos(2 * x); }, 0, M_PI, 1e-6);

    for (auto &interval: std::vector<std::pair<int, int>>
            {{1,  8},
             {9,  14},
             {15, 20},
             {29, 40}}) {
        int low, high;
        tie(low, high) = interval;

        int i, m;
        double E;
        auto eigenvalues = matslise.eigenvaluesByIndex(low, high);
        for (auto eigenvalue: eigenvalues) {
            tie(i, E, m) = eigenvalue;
            for (int j = i; j < i + m; ++j) {
                REQUIRE(seen.insert(j).second);
                auto f = exact.find(j);
                if (f != exact.end()) {
                    REQUIRE(Approx(f->second).epsilon(1e-4) == E);
                }
            }
        }
    }

    for (auto &ie: exact)
        REQUIRE(seen.find(ie.first) != seen.end());
}
