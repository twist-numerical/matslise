#include <cmath>
#include <set>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "../../src/util/quadrature.h"
#include "checkOrthonormality.h"
#include "../../src/util/constants.h"


using namespace matslise;
using namespace std;
using namespace Eigen;


using namespace Catch::Matchers;

void compareEigenfunctions(
        const Matslise2D<> &p, double E, const vector<function<double(double, double)>> &exact) {
    int n = 50, m = 60;
    ArrayXd x = ArrayXd::LinSpaced(n, p.domain.getMin(0), p.domain.getMax(0));
    ArrayXd y = ArrayXd::LinSpaced(m, p.domain.getMin(1), p.domain.getMax(1));
    const std::vector<Eigenfunction2D<>> fs = p.eigenfunction(E);

    REQUIRE(exact.size() == fs.size());
    for (const Eigenfunction2D<> &f :fs) {
        ArrayXXd fxy = f(x, y);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) {
                if (i % 10 == 0 && j % 10 == 0)
                    REQUIRE(Approx(f(x[i], y[j])).margin(1e-10) == fxy(i, j));
                bool valid = false;
                for (const function<double(double, double)> &e : exact)
                    if (abs(abs(e(x[i], y[j])) - abs(fxy(i, j))) < 1e-5)
                        valid = true;
                CHECKED_ELSE(valid) {
                    std::stringstream ss;
                    ss << E << ": " << x[i] << ", " << y[j] << " : " << f(x[i], x[j]);
                    FAIL(ss.str());
                }
            }
    }
}

TEST_CASE("Eigenvalues V=0", "[matslise2d][eigenfunctions][zero]") {
    Matslise2D<> p(
            [](double, double) -> double {
                return 0;
            },
            {{0, constants<double>::PI}, 0, constants<double>::PI},
            Options2<>().tolerance(1e-6).N(12).nested(Options1<>().tolerance(1e-5)));

    set<double> eigenvalues;
    vector<double> eigenvaluesList;
    for (int i = 1; i < 6; ++i) {
        for (int j = 1; j <= 6; ++j) {
            double E = i * i + j * j;
            if (eigenvalues.find(E) != eigenvalues.end())
                continue;
            eigenvaluesList.emplace_back(E);
            eigenvalues.insert(E);
            for (int k = -1; k <= 1; ++k)
                CHECK(Approx(p.eigenvalue(E + k * 1e-2).first).margin(1e-7) == E);

            vector<function<double(double, double)>> v;
            for (int k = 1; k * k < E; ++k) {
                int l = (int) round(sqrt(E - k * k));
                if (l * l == E - k * k) {
                    v.emplace_back([k, l](double x, double y) -> double {
                        return 2 * sin(x * k) * sin(y * l) / constants<double>::PI;
                    });
                }
            }
            compareEigenfunctions(p, E, v);
        }
    }
    sort(eigenvaluesList.begin(), eigenvaluesList.end());

    int i = 0;
    for (auto &iEm : p.eigenvaluesByIndex(0, 10)) {
        CHECK(Approx(get<1>(iEm)).margin(1e-7) == eigenvaluesList[i]);
        ++i;
    }

    vector<double> eigenvalues_v(eigenvalues.begin(), eigenvalues.end());
    checkOrthonormality(&p, eigenvalues_v.begin(), eigenvalues_v.end());
}
