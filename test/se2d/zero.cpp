#include <cmath>
#include <set>
#include <tuple>
#include "../catch.hpp"
#include "../../src/matslise.h"
#include "../../src/util/lobatto.h"
#include "checkOrthonormality.h"


using namespace matslise;
using namespace std;
using namespace Eigen;


using namespace Catch::Matchers;

void compareEigenfunctions(
        const SE2D<> &p, double E, const vector<function<double(double, double)>> &exact) {
    int n = 50, m = 60;
    ArrayXd x = ArrayXd::LinSpaced(n, p.domain.getMin(0), p.domain.getMax(0));
    ArrayXd y = ArrayXd::LinSpaced(m, p.domain.getMin(1), p.domain.getMax(1));
    const std::vector<ArrayXXd> fs = p.eigenfunction(E, x, y);

    REQUIRE(exact.size() == fs.size());
    for (const ArrayXXd &f :fs) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) {
                bool valid = false;
                for (const function<double(double, double)> &e : exact)
                    if (abs(abs(e(x[i], y[j])) - abs(f(i, j))) < 1e-5)
                        valid = true;
                CHECKED_ELSE(valid) {
                    std::stringstream ss;
                    ss << E << ": " << x[i] << ", " << y[j] << " : " << f(i, j);
                    FAIL(ss.str());
                }
            }
    }
}

TEST_CASE("Eigenvalues V=0", "[se2d][eigenfunctions][zero]") {
    SE2D<> p(
            [](double, double) -> double {
                return 0;
            },
            {{0, constants<double>::PI}, 0, constants<double>::PI},
            Options2<>().tolerance(1e-5).stepsPerSector(4).N(12).nested(Options1<>().tolerance(1e-5)));

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
                CHECK(Approx(p.eigenvalue(E + k * 1e-2)).margin(1e-7) == E);

            vector<function<double(double, double)>> v;
            for (int k = 1; k * k < E; ++k) {
                int l = (int) round(sqrt(E - k * k));
                if (l * l == E - k * k) {
                    v.emplace_back([k, l](double x, double y) -> double {
                        return 2 * sin(x * k) * sin(y * l) / constants<double>::PI;
                    });
                }
            }
            // can be executed when normalizing is done
            compareEigenfunctions(p, E, v);
        }
    }
    sort(eigenvaluesList.begin(), eigenvaluesList.end());

    int i = 0;
    for (const double E : p.firstEigenvalues(10)) {
        CHECK(Approx(E).margin(1e-7) == eigenvaluesList[i]);
        ++i;
    }

    vector<double> eigenvalues_v(eigenvalues.begin(), eigenvalues.end());
    checkOrthonormality(p, eigenvalues_v.begin(), eigenvalues_v.end());
}
