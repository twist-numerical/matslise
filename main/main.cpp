#include <iostream>
#include <cmath>
#include <set>
#include "../src/matslise.h"
#include "../src/matscs.h"
#include "../src/se2d.h"


#include <boost/multiprecision/float128.hpp>

using boost::multiprecision::float128;

using namespace std;
using namespace Eigen;
using namespace matslise;

void coffey() {
    double M = M_PI_2;
    double B = 20;

    matslise::Matslise<> coffey([B](double x) {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, -M, M, 128);

    cout << endl;

    matslise::Y<> y0 = Y<>::Dirichlet();
    ArrayXd x(6);
    x << 0.0700000000000000,
            0.0760000000000000,
            0.0820000000000000,
            0.0880000000000000,
            0.0940000000000000,
            0.100000000000000;

    vector<pair<int, double>> *eigenvalues = coffey.computeEigenvaluesByIndex(0, 2, y0, y0);
    cout.precision(17);
    unsigned int i;
    double E;
    for (pair<int, double> iE : *eigenvalues) {
        tie(i, E) = iE;

        Array<matslise::Y<>, Dynamic, 1> y = coffey.computeEigenfunction(E, y0, y0, x);
        for (int j = 0; j < x.size(); ++j) {
            cout << x[j] << ", ";
        }
        cout << endl;
        for (int j = 0; j < y.size(); ++j) {
            cout << y[j].y[0] << ", ";
        }

        cout << endl << i << ": " << E << endl << endl;
    }
    delete eigenvalues;
}

void mathieu() {
    double m = 0, M = M_PI;

    Matslise<> mathieu([](double x) {
        return 2 * cos(2 * x);
        //return -2 * 30 * cos(2 * x) + 30 * 30 * sin(2 * x) * sin(2 * x);
    }, m, M, Matslise<>::UNIFORM(16));
    double h = M_PI / 16;

    matslise::Y<> y0({0, 1}, {0, 0});
    cout << get<0>(mathieu.propagate(110, y0, m, m + 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(120, y0, m, m + 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(130, y0, m, m + 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(140, y0, m, m + 2 * h)) << endl;
    cout << endl;
    cout << get<0>(mathieu.propagate(110, y0, M, M - 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(120, y0, M, M - 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(130, y0, M, M - 2 * h)) << endl;
    cout << get<0>(mathieu.propagate(140, y0, M, M - 2 * h)) << endl;
    cout << endl;

    vector<pair<int, double>> *eigenvalues = mathieu.computeEigenvalues(0, 100, y0, y0);
    cout.precision(17);

    unsigned int i;
    double E;
    for (auto &iE : *eigenvalues) {
        tie(i, E) = iE;
        cout << i << ": " << E << endl;
    }
    delete eigenvalues;
}


void test2d() {
    std::cout.precision(10);
    SEnD<2> se2d([](double x, double y) { return (1 + x * x) * (1 + y * y); }, {{-5.5, 5.5}, -5.5, 5.5},
                 Options<2>().tolerance(1e-3).nested(Options<1>().sectorCount(26)).N(12));
    cout << "sectorCount: " << se2d.sectorCount << endl;
    cout << "match: " << se2d.match << endl;
    for (int i = 0; i < se2d.sectorCount; ++i)
        cout << "minmax: (" << se2d.sectors[i]->min << ", " << se2d.sectors[i]->max << "), vbar: "
             << se2d.sectors[i]->vbar.minCoeff() << endl;
    double E = se2d.findEigenvalue(8);
    cout << E << endl;

    ArrayXd xs(5);
    xs << -2, -1, 0, 1, 2;
    ArrayXd ys(5);
    ys << -2, -1, 0, 1, 2;
    vector<ArrayXXd> *fs = se2d.computeEigenfunction(E, {xs, ys});
    for (auto &a : *fs)
        cout << a << "\n\n" << endl;
    delete fs;
}

void testHenon() {
    SEnD<2> p(
            [](double x, double y) -> double {
                return (x * x + y * y) + 1 / (2 * sqrt(5)) * x * (y * y - x * x / 3);
            },
            {{-6, 6}, -6, 6},
            Options<2>().stepsPerSector(1).sectorCount(13).nested(Options<1>().sectorCount(16)));

    pair<double, int> eigenvalues[] = {
            {0.998594690530479, 1},
            {1.99007660445524,  2},
            {2.95624333869018,  1},
            {2.98532593386986,  2},
            {3.92596412795287,  2},
            {3.98241882458866,  1},
            {3.98575763690663,  1},
            {4.87014557482289,  1},
            {4.89864497284387,  2}};

    ArrayXd x(3);
    x << -1, 0, 1;
    double E;
    int multiplicity;
    for (auto &Em  : eigenvalues) {
        tie(E, multiplicity) = Em;
        E *= 2;

        const pair<double, double> &error = p.calculateError(E);
        cout << error.first << ", " << error.second << endl;

        const vector<Array<double, -1, -1>> *f = p.computeEigenfunction(E, {x, x});
        cout << f->size() << " == " << multiplicity << endl;
        delete f;
    }
}


void compareEigenfunctions(
        const SEnD<2> &p, double E, const vector<function<double(double, double)>> &exact) {
    int n = 50, m = 60;
    ArrayXd x = ArrayXd::LinSpaced(n, p.domain.getMin(0), p.domain.getMax(0));
    ArrayXd y = ArrayXd::LinSpaced(m, p.domain.getMin(1), p.domain.getMax(1));
    std::vector<ArrayXXd> *fs = p.computeEigenfunction(E, {x, y});

    cout << exact.size() << " == " << fs->size();
    for (ArrayXXd &f :*fs) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) {
                bool valid = false;
                for (const function<double(double, double)> &e : exact)
                    if (abs(abs(e(x[i], y[j])) - abs(f(i, j))) < 1e-7)
                        valid = true;
                if (!valid) {
                    cout << E << ": " << x[i] << ", " << y[j] << " : " << f(i, j) << endl;
                }
            }
    }
    delete fs;
}

void testZero() {
    SEnD<2> p(
            [](double x, double y) -> double {
                return 0;
            },
            {{0, M_PI}, 0, M_PI},
            Options<2>().sectorCount(13).stepsPerSector(4).N(12).nested(Options<1>().sectorCount(13)));

    set<double> eigenvalues;
    for (int i = 1; i < 6; ++i) {
        for (int j = 1; j <= i; ++j) {
            double E = i * i + j * j;
            if (eigenvalues.find(E) != eigenvalues.end())
                continue;
            eigenvalues.insert(E);
            for (int k = -1; k <= 1; ++k)
                cout << p.findEigenvalue(E + k * 1e-2) << "==" << E << endl;

            vector<function<double(double, double)>> v;
            for (int k = 1; k * k < E; ++k) {
                int l = (int) round(sqrt(E - k * k));
                if (l * l == E - k * k) {
                    v.push_back([k, l](double x, double y) -> double {
                        return sin(x * k) * sin(y * l) / M_PI_2;
                    });
                }
            }
            compareEigenfunctions(p, E, v);
        }
    }
}

/*
void test3d() {
    std::cout.precision(10);
    SEnD<3> se3d([](double x, double y, double z) -> double { return (1 + x * x) * (1 + y * y) * (1 + z * z); },
                 {{{-5.5, 5.5}, -5.5, 5.5}, -5.5, 5.5},
                 Options<3>().sectorCount(10).stepsPerSector(5).N(8).nested(
                         Options<2>().sectorCount(25).stepsPerSector(5).N(8).nested(
                                 Options<1>().sectorCount(35))));

    for (double a = 3.1; a < 3.3; a += .01)
        cout << se3d.calculateError(a).first << endl;
}*/

void testMatscs() {
    Matscs<> ms([](double x) {
        return (MatrixXd(1, 1) << x * x).finished();
    }, 1, -5, 5, 32);
    MatrixXd zero = MatrixXd::Zero(1, 1);
    MatrixXd one = MatrixXd::Identity(1, 1);
    double E = 3;
    cout << ms.propagate(E, Y<double, Dynamic>::Dirichlet(1), -5, 5).y(0, 0) << endl;

    cout << "\n\nMATSLISE" << endl;
    Matslise<> ml([](double x) { return x * x; }, -5, 5, 32);
    cout << get<0>(ml.propagate(E, Y<>({0, 1}, {0, 0}), -5, 5)).y(0, 0) << endl;
}

void testEigenfunctionCalculator() {
    double M = M_PI_2;
    double B = 20;

    Matslise<> coffey([B](double x) {
        return -2 * B * cos(2 * x) + B * B * sin(2 * x) * sin(2 * x);
    }, -M, M, 128);

    Y<> y0({0, 1}, {0, 0});

    vector<pair<int, double>> *eigs = coffey.computeEigenvaluesByIndex(2, 3, y0, y0);
    double E = get<1>((*eigs)[0]);
    cout << E << endl;
    delete eigs;

    function<Y<>(double)> ef = coffey.eigenfunctionCalculator(E, y0, y0);
    cout << ef(0) << endl;
}

void testPrufer() {
    //Matslise ms([](double x) { return 0.0; }, 0., 3.14, 60);

    Matslise<> ms([](double x) {
        return 2 * 5 * cos(2 * x);
    }, -M_PI_2, M_PI_2, 18);
    Y<> y0({0, 1});
    cout << get<1>(ms.propagate(1, y0, -M_PI_2, -1)) << endl;
    // cout << get<0>(ms.propagate(16, y0, .1, 1.6)) << endl;
}

void testHigh() {

    Matslise<> ms([](double x) -> double {
        return (1 - cos(2 * M_PI * x)) / 2 * 1000;
    }, 0, 1, Matslise<>::AUTO(1e-6));
    cout << ms.computeEigenvaluesByIndex(0, 10, Y<>::Dirichlet(), Y<>::Dirichlet()) << endl;
}

void testQuartic() {
    float128 a = -1q;
    float128 c = 0.001q;
    SEnD<2>([a, c](double x, double y) -> double  {
        return x*x+y*y+c*(x*x*x*x+2*a*x*x*y*y+y*y*y*y);
    }, {-7,7});
}

int main() {
    // coffey();
    // testPrufer();
    // test2d();
    // testHenon();
    // testZero();
    // test3d();
    // testBigE();
    // testMatscs();
    // testEigenfunctionCalculator();
    // testHigh();
    mathieu();
}