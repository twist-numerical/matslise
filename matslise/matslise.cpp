//
// Created by toon on 5/16/18.
//

#include <cmath>
#include <array>
#include <vector>
#include <queue>
#include "matslise.h"
#include "legendre.h"
#include "calculateEta.h"

#define EPS (1.e-12)

using namespace matslise;
using namespace std;

Matslise::Matslise(std::function<double(double)> V, double xmin, double xmax, int sectorCount)
        : V(V), xmin(xmin), xmax(xmax), sectorCount(sectorCount) {
    sectors = new Sector *[sectorCount];
    double h = (xmax - xmin) / sectorCount;
    double mid = (xmax + xmin) / 2;
    for (int i = 0; i < sectorCount; ++i) {
        double a = xmin + i * h;
        double b = a + h;
        if (b - 1.e-5 > mid) {
            match = b;
            mid = xmax + 1;
        }
        sectors[i] = new Sector(this, a, b);
    }
}


std::tuple<Y, double> Matslise::propagate(double E, const Y &_y, double a, double b) const {
    Y y = _y;
    double theta = y.theta();
    if (a < b) {
        for (int i = 0; i < sectorCount; ++i) {
            Sector *sector = sectors[i];
            if (sector->xmax > a) {
                if (sector->xmin < a) // first
                    y = sector->propagate(E, y, sector->xmin - a, theta);

                if (sector->xmax > b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, theta);
                    break;
                }

                y = sector->propagate(E, y, sector->xmax - sector->xmin, theta);
            }
        }
    } else {
        for (int i = sectorCount - 1; i >= 0; --i) {
            Sector *sector = sectors[i];
            if (sector->xmin < a) {
                if (sector->xmax > a) // first
                    y = sector->propagate(E, y, sector->xmin - a, theta);
                else
                    y = sector->propagate(E, y, sector->xmin - sector->xmax, theta);

                if (sector->xmin < b) { // last
                    y = sector->propagate(E, y, b - sector->xmin, theta);
                    break;
                }

            }
        }
    }
    return make_tuple(y, theta);
}

tuple<double, double, double>
Matslise::calculateError(double E, const Y &left, const Y &right) const {
    Y l, r;
    double thetaL, thetaR;
    tie(l, thetaL) = propagate(E, left, xmin, match);
    tie(r, thetaR) = propagate(E, right, xmax, match);
    return make_tuple(l.y[1] * r.y[0] - r.y[1] * l.y[0],
                           l.dy[1] * r.y[0] + l.y[1] * r.dy[0] - (r.dy[1] * l.y[0] + r.y[1] * l.dy[0]),
                           thetaL - thetaR);
}

double newtonIteration(const Matslise *ms, double E, const Y &left, const Y &right, double tol) {
    double adjust, error, derror, theta;
    int i = 0;
    do {
        tie(error, derror, theta) = ms->calculateError(E, left, right);
        adjust = error/derror;
        E = E - adjust;
        if(++i > 50) {
            throw runtime_error("Newton-iteration did not converge");
        }
    } while(fabs(adjust) > tol);
    return E;
}

vector<double> *Matslise::computeEigenvalues(double Emin, double Emax, const Y &left, const Y &right) const {
    vector<double> *eigenvalues =  new vector<double>();
    queue<tuple<double, double, double, double>> toCheck;

    toCheck.push(make_tuple(Emin, get<2>(calculateError(Emin, left, right))/M_PI,
            Emax, get<2>(calculateError(Emax, left, right))/M_PI));

    double a, ta, b, tb, c, tc;
    int ia, ib;
    while(!toCheck.empty()) {
        tie(a, ta, b, tb) = toCheck.front();
        toCheck.pop();
        ia = (int) floor(ta);
        ib = (int) floor(tb);
        if(ta >= tb || ia == ib)
            continue;

        c = (a+b)/2;
        if(tb - ta < 0.1)
            eigenvalues->push_back(newtonIteration(this, c, left, right, 1e-12));
        else {
            tc = get<2>(calculateError(c, left, right))/M_PI;
            toCheck.push(make_tuple(a, ta, c, tc));
            toCheck.push(make_tuple(c, tc, b, tb));
        }
    }

    sort(eigenvalues->begin(), eigenvalues->end());

    return eigenvalues;
}

Matslise::~Matslise() {
    for (int i = 0; i < sectorCount; ++i)
        delete sectors[i];
    delete[] sectors;
}

std::vector<Y> *Matslise::computeEigenfunction(double E, std::vector<double> &x) const {
    std::sort(x.begin(), x.end());
    std::vector<Y> *ys = new std::vector<Y>();

    auto iterator = x.begin();

    while (iterator != x.end() && *iterator < xmin - EPS)
        iterator = x.erase(iterator);

    Sector *sector;
    Y y({0, 1});
    for (int i = 0; iterator != x.end(); ++iterator) {
        while ((sector = sectors[i])->xmax < *iterator) {
            y = sector->calculateT(E) * y;
            ++i;
            if (i >= sectorCount)
                goto allSectorsDone;
        }

        ys->push_back(sector->calculateT(E, *iterator - sector->xmin) * y);
    }
    allSectorsDone:

    while (iterator != x.end() && *iterator > xmax + EPS)
        iterator = x.erase(iterator);

    return ys;
}

Sector::Sector(Matslise *s, double xmin, double xmax) : s(s), xmin(xmin), xmax(xmax) {
    h = xmax - xmin;
    vs = legendre::getCoefficients(6, s->V, xmin, xmax);

    calculateTCoeffs();
}

void Sector::calculateTCoeffs() {
    double v1 = vs[1],
            v2 = vs[2],
            v3 = vs[3],
            v4 = vs[4],
            v5 = vs[5];

    // @formatter:off
    u = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0, 0.0, (((((((h*h*h*h*h)*v5)*-0.5)+(((h*h*h*h)*v4)*0.5))+(((h*h*h)*v3)*-0.5))+(((h*h)*v2)*0.5))+((h*v1)*-0.5)), (((((((h*h*h*h)*v5)*7.5)+(((h*h*h)*v4)*-5.0))+(((h*h)*v3)*3.0))+((h*v2)*-1.5))+(v1*0.5)), ((((((h*h*h)*v5)*-35.0)+(((h*h)*v4)*15.0))+((h*v3)*-5.0))+v2), (((((h*h)*v5)*70.0)+((h*v4)*-17.5))+(v3*2.5)), (((h*v5)*-63.0)+(v4*7.0)), (v5*21.0),   0.0, 0.0, 0.0, (((((((h*h*h*h)*v5)*-7.5)+(((h*h*h)*v4)*5.0))+(((h*h)*v3)*-3.0))+((h*v2)*1.5))+(v1*-0.5)), (((((((((((v3*v4)+(v2*v5))*(h*h*h*h*h*h*h))*-0.25)+(((((v3*v3)+((v2*v4)*2.0))+((v1*v5)*2.0))*(h*h*h*h*h*h))*0.125))+((((v2*v3)+(v1*v4))*(h*h*h*h*h))*-0.25))+((((v2*v2)+((v1*v3)*2.0))*(h*h*h*h))*0.125))+((((v1*v2)+(v5*-210.0))*(h*h*h))*-0.25))+((((v1*v1)+(v4*-180.0))*(h*h))*0.125))+((h*v3)*7.5))+(v2*-1.5)), (((((((((((v3*v4)*8.0)+((v2*v5)*9.0))*(h*h*h*h*h*h))*0.5)+((((((v3*v3)*6.0)+((v2*v4)*13.0))+((v1*v5)*16.0))*(h*h*h*h*h))*-0.25))+(((((v2*v3)*9.0)+((v1*v4)*11.0))*(h*h*h*h))*0.25))+(((((v2*v2)*3.0)+((v1*v3)*7.0))*(h*h*h))*-0.25))+(((v1*v2)+(v5*-140.0))*(h*h)))+((((v1*v1)+(v4*-140.0))*h)*-0.25))+(v3*-5.0)), (((((((((((v3*v4)*100.0)+((v2*v5)*117.0))*(h*h*h*h*h))*-0.25)+((((((v3*v3)*28.0)+((v2*v4)*62.0))+((v1*v5)*85.0))*(h*h*h*h))*0.25))+(((((v2*v3)*3.0)+((v1*v4)*4.0))*(h*h*h))*-2.5))+(((((v2*v2)*13.0)+((v1*v3)*32.0))*(h*h))*0.125))+((((v1*v2)+(v5*-126.0))*h)*-1.25))+((v1*v1)*0.125))+(v4*-17.5)), ((((((((((v3*v4)*16.0)+((v2*v5)*19.0))*(h*h*h*h))*5.0)+((((((v3*v3)*13.0)+((v2*v4)*29.0))+((v1*v5)*42.0))*(h*h*h))*-1.25))+(((((v2*v3)*47.0)+((v1*v4)*65.0))*(h*h))*0.25))+(((((v2*v2)*2.0)+((v1*v3)*5.0))*h)*-0.75))+((v1*v2)*0.5))+(v5*-63.0)),   0.0, 0.0, 0.0, 0.0, 0.0, (((((((((((v3*v4)*8.0)+((v2*v5)*9.0))*(h*h*h*h*h*h))*-0.5)+((((((v3*v3)*6.0)+((v2*v4)*13.0))+((v1*v5)*16.0))*(h*h*h*h*h))*0.25))+(((((v2*v3)*9.0)+((v1*v4)*11.0))*(h*h*h*h))*-0.25))+(((((v2*v2)*3.0)+((v1*v3)*7.0))*(h*h*h))*0.25))+((((v1*v2)+(v5*-70.0))*(h*h))*-1.0))+((((v1*v1)+(v4*-70.0))*h)*0.25))+(v3*2.5)), ((((((((((((v1*(v2*v2))+((v1*v1)*v3))+((v3*v4)*-800.0))+((v2*v5)*-852.0))*(h*h*h*h*h))*-0.0625)+(((((((v1*v1)*v2)+((v3*v3)*-228.0))+((v2*v4)*-472.0))+((v1*v5)*-560.0))*(h*h*h*h))*0.0625))+(((((v1*v1*v1)+((v2*v3)*-720.0))+((v1*v4)*-820.0))*(h*h*h))*-0.0208333333333))+(((((v2*v2)*27.0)+((v1*v3)*58.0))*(h*h))*-0.125))+((((v1*v2)+(v5*-63.0))*h)*2.5))+((v1*v1)*-0.291666666667))+(v4*17.5)), ((((((((((((v1*(v2*v2))*7.0)+(((v1*v1)*v3)*8.0))+((v3*v4)*-3680.0))+((v2*v5)*-4000.0))*(h*h*h*h))*0.0625)+(((((((v1*v1)*v2)+((v3*v3)*-152.0))+((v2*v4)*-320.0))+((v1*v5)*-392.0))*(h*h*h))*-0.3125))+(((((v1*v1*v1)+((v2*v3)*-544.0))+((v1*v4)*-640.0))*(h*h))*0.0625))+(((((v2*v2)*9.0)+((v1*v3)*20.0))*h)*0.5))+((v1*v2)*-1.5))+(v5*94.5)),   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ((((((((((((v1*(v2*v2))*7.0)+(((v1*v1)*v3)*8.0))+((v3*v4)*-1280.0))+((v2*v5)*-1520.0))*(h*h*h*h))*-0.0625)+(((((((v1*v1)*v2)+((v3*v3)*-52.0))+((v2*v4)*-116.0))+((v1*v5)*-168.0))*(h*h*h))*0.3125))+(((((v1*v1*v1)+((v2*v3)*-188.0))+((v1*v4)*-260.0))*(h*h))*-0.0625))+(((((v2*v2)*2.0)+((v1*v3)*5.0))*h)*-0.75))+((v1*v2)*0.5))+(v5*-31.5))};
    up = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0, (((((((h*h*h*h*h)*v5)*-0.5)+(((h*h*h*h)*v4)*0.5))+(((h*h*h)*v3)*-0.5))+(((h*h)*v2)*0.5))+((h*v1)*-0.5)), (((((((h*h*h*h)*v5)*7.5)+(((h*h*h)*v4)*-5.0))+(((h*h)*v3)*3.0))+((h*v2)*-1.5))+(v1*0.5)), (((((((((((v3*v4)+(v2*v5))*(h*h*h*h*h*h*h))*-0.25)+(((((v3*v3)+((v2*v4)*2.0))+((v1*v5)*2.0))*(h*h*h*h*h*h))*0.125))+((((v2*v3)+(v1*v4))*(h*h*h*h*h))*-0.25))+((((v2*v2)+((v1*v3)*2.0))*(h*h*h*h))*0.125))+((((v1*v2)+(v5*210.0))*(h*h*h))*-0.25))+((((v1*v1)+(v4*180.0))*(h*h))*0.125))+((h*v3)*-7.5))+(v2*1.5)), (((((((((((v3*v4)*8.0)+((v2*v5)*9.0))*(h*h*h*h*h*h))*0.5)+((((((v3*v3)*6.0)+((v2*v4)*13.0))+((v1*v5)*16.0))*(h*h*h*h*h))*-0.25))+(((((v2*v3)*9.0)+((v1*v4)*11.0))*(h*h*h*h))*0.25))+(((((v2*v2)*3.0)+((v1*v3)*7.0))*(h*h*h))*-0.25))+(((v1*v2)+(v5*140.0))*(h*h)))+((((v1*v1)+(v4*140.0))*h)*-0.25))+(v3*5.0)), (((((((((((v3*v4)*100.0)+((v2*v5)*117.0))*(h*h*h*h*h))*-0.25)+((((((v3*v3)*28.0)+((v2*v4)*62.0))+((v1*v5)*85.0))*(h*h*h*h))*0.25))+(((((v2*v3)*3.0)+((v1*v4)*4.0))*(h*h*h))*-2.5))+(((((v2*v2)*13.0)+((v1*v3)*32.0))*(h*h))*0.125))+((((v1*v2)+(v5*126.0))*h)*-1.25))+((v1*v1)*0.125))+(v4*17.5)), ((((((((((v3*v4)*16.0)+((v2*v5)*19.0))*(h*h*h*h))*5.0)+((((((v3*v3)*13.0)+((v2*v4)*29.0))+((v1*v5)*42.0))*(h*h*h))*-1.25))+(((((v2*v3)*47.0)+((v1*v4)*65.0))*(h*h))*0.25))+(((((v2*v2)*2.0)+((v1*v3)*5.0))*h)*-0.75))+((v1*v2)*0.5))+(v5*63.0)), (((((((((v3*v4)*41.0)+((v2*v5)*49.0))*(h*h*h))*-3.5)+((((((v3*v3)*80.0)+((v2*v4)*179.0))+((v1*v5)*266.0))*(h*h))*0.25))+(((((v2*v3)*5.0)+((v1*v4)*7.0))*h)*-1.75))+((v2*v2)*0.5))+((v1*v3)*1.25)),   0.0, 0.0, 0.0, (((((((((((v3*v4)+(v2*v5))*(h*h*h*h*h*h*h))*-0.25)+(((((v3*v3)+((v2*v4)*2.0))+((v1*v5)*2.0))*(h*h*h*h*h*h))*0.125))+((((v2*v3)+(v1*v4))*(h*h*h*h*h))*-0.25))+((((v2*v2)+((v1*v3)*2.0))*(h*h*h*h))*0.125))+((((v1*v2)+(v5*-210.0))*(h*h*h))*-0.25))+((((v1*v1)+(v4*-180.0))*(h*h))*0.125))+((h*v3)*7.5))+(v2*-1.5)), (((((((((((v3*v4)*8.0)+((v2*v5)*9.0))*(h*h*h*h*h*h))*0.5)+((((((v3*v3)*6.0)+((v2*v4)*13.0))+((v1*v5)*16.0))*(h*h*h*h*h))*-0.25))+(((((v2*v3)*9.0)+((v1*v4)*11.0))*(h*h*h*h))*0.25))+(((((v2*v2)*3.0)+((v1*v3)*7.0))*(h*h*h))*-0.25))+(((v1*v2)+(v5*-210.0))*(h*h)))+((((v1*v1)+(v4*-210.0))*h)*-0.25))+(v3*-7.5)), ((((((((((((v1*(v2*v2))+((v1*v1)*v3))+((v3*v4)*400.0))+((v2*v5)*552.0))*(h*h*h*h*h))*-0.0625)+(((((((v1*v1)*v2)+((v3*v3)*108.0))+((v2*v4)*272.0))+((v1*v5)*460.0))*(h*h*h*h))*0.0625))+(((((v1*v1*v1)+((v2*v3)*360.0))+((v1*v4)*620.0))*(h*h*h))*-0.0208333333333))+(((((v2*v2)*6.0)+((v1*v3)*19.0))*(h*h))*0.25))+((((v1*v2)+(v5*-252.0))*h)*-1.25))+((v1*v1)*0.0833333333333))+(v4*-35.0)), ((((((((((((v1*(v2*v2))*7.0)+(((v1*v1)*v3)*8.0))+((v3*v4)*1440.0))+((v2*v5)*2080.0))*(h*h*h*h))*0.0625)+(((((((v1*v1)*v2)+((v3*v3)*56.0))+((v2*v4)*144.0))+((v1*v5)*280.0))*(h*h*h))*-0.3125))+(((((v1*v1*v1)+((v2*v3)*208.0))+((v1*v4)*400.0))*(h*h))*0.0625))+(((((v2*v2)*3.0)+((v1*v3)*10.0))*h)*-0.5))+((v1*v2)*0.5))+(v5*-157.5)), (((((((((((v1*(v2*v2))*19.0)+(((v1*v1)*v3)*23.0))+((v3*v4)*3036.0))+((v2*v5)*4284.0))*(h*h*h))*-0.0625)+((((((((v1*v1)*v2)*9.0)+((v3*v3)*408.0))+((v2*v4)*1016.0))+((v1*v5)*2044.0))*(h*h))*0.0625))+(((((v1*v1*v1)+((v2*v3)*180.0))+((v1*v4)*336.0))*h)*-0.0625))+((v2*v2)*0.6))+((v1*v3)*1.75)),   0.0, 0.0, 0.0, 0.0, 0.0, ((((((((((((v1*(v2*v2))+((v1*v1)*v3))+((v3*v4)*-800.0))+((v2*v5)*-852.0))*(h*h*h*h*h))*-0.0625)+(((((((v1*v1)*v2)+((v3*v3)*-228.0))+((v2*v4)*-472.0))+((v1*v5)*-560.0))*(h*h*h*h))*0.0625))+(((((v1*v1*v1)+((v2*v3)*-720.0))+((v1*v4)*-820.0))*(h*h*h))*-0.0208333333333))+(((((v2*v2)*27.0)+((v1*v3)*58.0))*(h*h))*-0.125))+((((v1*v2)+(v5*-63.0))*h)*2.5))+((v1*v1)*-0.291666666667))+(v4*17.5)), ((((((((((((v1*(v2*v2))*7.0)+(((v1*v1)*v3)*8.0))+((v3*v4)*-6080.0))+((v2*v5)*-6480.0))*(h*h*h*h))*0.0625)+(((((((v1*v1)*v2)+((v3*v3)*-252.0))+((v2*v4)*-524.0))+((v1*v5)*-616.0))*(h*h*h))*-0.3125))+(((((v1*v1*v1)+((v2*v3)*-900.0))+((v1*v4)*-1020.0))*(h*h))*0.0625))+(((((v2*v2)*6.0)+((v1*v3)*13.0))*h)*1.25))+((v1*v2)*-2.5))+(v5*157.5)), (((((((((((v1*(v2*v2))*48.0)+(((v1*v1)*v3)*71.0))+((v3*v4)*-55332.0))+((v2*v5)*-60228.0))*(h*h*h))*-0.0208333333333)+((((((((v1*v1)*v2)*23.0)+((v3*v3)*-7866.0))+((v2*v4)*-16572.0))+((v1*v5)*-20328.0))*(h*h))*0.0208333333333))+(((((v1*v1*v1)+((v2*v3)*-1710.0))+((v1*v4)*-2016.0))*h)*-0.0416666666667))+((v2*v2)*-4.275))+((v1*v3)*-9.5)),   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (((((((((((v1*(v2*v2))*123.0)+(((v1*v1)*v3)*136.0))+((v3*v4)*-20664.0))+((v2*v5)*-24696.0))*(h*h*h))*0.0208333333333)+((((((((v1*v1)*v2)*29.0)+((v3*v3)*-1431.0))+((v2*v4)*-3252.0))+((v1*v5)*-4578.0))*(h*h))*-0.0416666666667))+(((((v1*v1*v1)+((v2*v3)*-180.0))+((v1*v4)*-252.0))*h)*0.145833333333))+((v2*v2)*1.425))+((v1*v3)*4.0))};
    v = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0, 0.0, 0.0, (((((((h*h*h*h*h)*v5)*-0.5)+(((h*h*h*h)*v4)*0.5))+(((h*h*h)*v3)*-0.5))+(((h*h)*v2)*0.5))+((h*v1)*-0.5)), (((((((h*h*h*h)*v5)*7.5)+(((h*h*h)*v4)*-5.0))+(((h*h)*v3)*3.0))+((h*v2)*-1.5))+(v1*0.5)), ((((((h*h*h)*v5)*-35.0)+(((h*h)*v4)*15.0))+((h*v3)*-5.0))+v2), (((((h*h)*v5)*70.0)+((h*v4)*-17.5))+(v3*2.5)), (((h*v5)*-63.0)+(v4*7.0)),   0.0, 0.0, 0.0, 0.0, 0.0, (((((((((((v3*v3)+((v2*v4)*2.0))+((v1*v5)*2.0))*(h*h*h*h*h*h))*0.125)+((((v2*v3)+(v1*v4))*(h*h*h*h*h))*-0.25))+((((v2*v2)+((v1*v3)*2.0))*(h*h*h*h))*0.125))+((((v1*v2)+(v5*-70.0))*(h*h*h))*-0.25))+((((v1*v1)+(v4*-60.0))*(h*h))*0.125))+((h*v3)*2.5))+(v2*-0.5)), (((((((((((v3*v3)*6.0)+((v2*v4)*13.0))+((v1*v5)*16.0))*(h*h*h*h*h))*-0.25)+(((((v2*v3)*9.0)+((v1*v4)*11.0))*(h*h*h*h))*0.25))+(((((v2*v2)*3.0)+((v1*v3)*7.0))*(h*h*h))*-0.25))+(((v1*v2)+(v5*-70.0))*(h*h)))+((((v1*v1)+(v4*-70.0))*h)*-0.25))+(v3*-2.5)), (((((((((((v3*v3)*28.0)+((v2*v4)*62.0))+((v1*v5)*85.0))*(h*h*h*h))*0.25)+(((((v2*v3)*3.0)+((v1*v4)*4.0))*(h*h*h))*-2.5))+(((((v2*v2)*13.0)+((v1*v3)*32.0))*(h*h))*0.125))+(((((v1*v2)*5.0)+(v5*-378.0))*h)*-0.25))+((v1*v1)*0.125))+(v4*-10.5)),   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ((((((((((((v1*v1)*v2)+((v3*v3)*-44.0))+((v2*v4)*-104.0))+((v1*v5)*-160.0))*(h*h*h*h))*0.0625)+(((((v1*v1*v1)+((v2*v3)*-144.0))+((v1*v4)*-220.0))*(h*h*h))*-0.0208333333333))+(((((v2*v2)*5.0)+((v1*v3)*14.0))*(h*h))*-0.125))+((((v1*v2)+(v5*-63.0))*h)*0.5))+((v1*v1)*-0.0416666666667))+(v4*3.5))};
    vp = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0, 0.0, (((((((h*h*h*h*h)*v5)*-0.5)+(((h*h*h*h)*v4)*0.5))+(((h*h*h)*v3)*-0.5))+(((h*h)*v2)*0.5))+((h*v1)*-0.5)), (((((((h*h*h*h)*v5)*7.5)+(((h*h*h)*v4)*-5.0))+(((h*h)*v3)*3.0))+((h*v2)*-1.5))+(v1*0.5)), ((((((h*h*h)*v5)*-35.0)+(((h*h)*v4)*15.0))+((h*v3)*-5.0))+v2), (((((h*h)*v5)*70.0)+((h*v4)*-17.5))+(v3*2.5)), (((h*v5)*-63.0)+(v4*7.0)), (v5*21.0),   0.0, 0.0, 0.0, (((((((h*h*h*h)*v5)*7.5)+(((h*h*h)*v4)*-5.0))+(((h*h)*v3)*3.0))+((h*v2)*-1.5))+(v1*0.5)), (((((((((((v3*v3)+((v2*v4)*2.0))+((v1*v5)*2.0))*(h*h*h*h*h*h))*0.125)+((((v2*v3)+(v1*v4))*(h*h*h*h*h))*-0.25))+((((v2*v2)+((v1*v3)*2.0))*(h*h*h*h))*0.125))+((((v1*v2)+(v5*210.0))*(h*h*h))*-0.25))+((((v1*v1)+(v4*180.0))*(h*h))*0.125))+((h*v3)*-7.5))+(v2*1.5)), (((((((((((v3*v3)*6.0)+((v2*v4)*13.0))+((v1*v5)*16.0))*(h*h*h*h*h))*-0.25)+(((((v2*v3)*9.0)+((v1*v4)*11.0))*(h*h*h*h))*0.25))+(((((v2*v2)*3.0)+((v1*v3)*7.0))*(h*h*h))*-0.25))+(((v1*v2)+(v5*140.0))*(h*h)))+((((v1*v1)+(v4*140.0))*h)*-0.25))+(v3*5.0)), (((((((((((v3*v3)*28.0)+((v2*v4)*62.0))+((v1*v5)*85.0))*(h*h*h*h))*0.25)+(((((v2*v3)*3.0)+((v1*v4)*4.0))*(h*h*h))*-2.5))+(((((v2*v2)*13.0)+((v1*v3)*32.0))*(h*h))*0.125))+((((v1*v2)+(v5*126.0))*h)*-1.25))+((v1*v1)*0.125))+(v4*17.5)), ((((((((((v3*v3)*13.0)+((v2*v4)*29.0))+((v1*v5)*42.0))*(h*h*h))*-1.25)+(((((v2*v3)*47.0)+((v1*v4)*65.0))*(h*h))*0.25))+(((((v2*v2)*2.0)+((v1*v3)*5.0))*h)*-0.75))+((v1*v2)*0.5))+(v5*63.0)),   0.0, 0.0, 0.0, 0.0, 0.0, (((((((((((v3*v3)*6.0)+((v2*v4)*13.0))+((v1*v5)*16.0))*(h*h*h*h*h))*-0.25)+(((((v2*v3)*9.0)+((v1*v4)*11.0))*(h*h*h*h))*0.25))+(((((v2*v2)*3.0)+((v1*v3)*7.0))*(h*h*h))*-0.25))+(((v1*v2)+(v5*-70.0))*(h*h)))+((((v1*v1)+(v4*-70.0))*h)*-0.25))+(v3*-2.5)), ((((((((((((v1*v1)*v2)+((v3*v3)*180.0))+((v2*v4)*392.0))+((v1*v5)*520.0))*(h*h*h*h))*0.0625)+(((((v1*v1*v1)+((v2*v3)*576.0))+((v1*v4)*740.0))*(h*h*h))*-0.0208333333333))+(((((v2*v2)*21.0)+((v1*v3)*50.0))*(h*h))*0.125))+(((((v1*v2)*4.0)+(v5*-315.0))*h)*-0.5))+((v1*v1)*0.208333333333))+(v4*-17.5)), (((((((((((v1*v1)*v2)+((v3*v3)*104.0))+((v2*v4)*232.0))+((v1*v5)*336.0))*(h*h*h))*-0.3125)+(((((v1*v1*v1)+((v2*v3)*376.0))+((v1*v4)*520.0))*(h*h))*0.0625))+(((((v2*v2)*2.0)+((v1*v3)*5.0))*h)*-1.5))+(v1*v2))+(v5*-94.5)),   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (((((((((((v1*v1)*v2)+((v3*v3)*-52.0))+((v2*v4)*-116.0))+((v1*v5)*-168.0))*(h*h*h))*-0.3125)+(((((v1*v1*v1)+((v2*v3)*-188.0))+((v1*v4)*-260.0))*(h*h))*0.0625))+(((((v2*v2)*2.0)+((v1*v3)*5.0))*h)*0.75))+((v1*v2)*-0.5))+(v5*31.5))};
    // @formatter:on

    for (int i = 0; i < MATSLISE_ETA; ++i) {
        hu[i] = hup[i] = hv[i] = hvp[i] = 0;
        double H = 1;
        for (int j = 0; j < MATSLISE_HMAX; ++j, H *= h) {
            hu[i] += H * u[i][j];
            hup[i] += H * up[i][j];
            hv[i] += H * v[i][j];
            hvp[i] += H * vp[i][j];
        }
    }
}

T Sector::calculateT(double E, double delta) const {
    if (fabs(delta) <= EPS)
        return T();
    if (fabs(delta - h) <= EPS)
        return calculateT(E);

    double *eta = calculateEta((vs[0] - E) * delta * delta, MATSLISE_ETA);
    T t((Matrix2d() << 0, 0, (vs[0] - E) * delta * eta[1], 0).finished(),
        (Matrix2d() << 0, 0, -delta * eta[1] + -(vs[0] - E) * delta * delta * delta * eta[2] / 2, 0).finished());

    for (int i = 0; i < MATSLISE_ETA; ++i) {
        double D = 1;
        for (int j = 0; j < MATSLISE_HMAX; ++j, D *= delta) {
            t.t(0, 0) += D * eta[i] * u[i][j];
            t.t(0, 1) += D * eta[i] * v[i][j];
            t.t(1, 0) += D * eta[i] * up[i][j];
            t.t(1, 1) += D * eta[i] * vp[i][j];

            if (i + 1 < MATSLISE_ETA) {
                double dEta = -delta * delta * eta[i + 1] / 2;
                t.dt(0, 0) += D * dEta * u[i][j];
                t.dt(0, 1) += D * dEta * v[i][j];
                t.dt(1, 0) += D * dEta * up[i][j];
                t.dt(1, 1) += D * dEta * vp[i][j];
            }
        }
    }

    delete[] eta;
    return t;
}

T Sector::calculateT(double E) const {
    double *eta = calculateEta((vs[0] - E) * h * h, MATSLISE_ETA);
    T t((Matrix2d() << 0, 0, (vs[0] - E) * h * eta[1], 0).finished(),
        (Matrix2d() << 0, 0, -h * eta[1] + -(vs[0] - E) * h * h * h * eta[2] / 2, 0).finished());

    for (int i = 0; i < MATSLISE_ETA; ++i) {
        t.t(0, 0) += eta[i] * hu[i];
        t.t(0, 1) += eta[i] * hv[i];
        t.t(1, 0) += eta[i] * hup[i];
        t.t(1, 1) += eta[i] * hvp[i];

        if (i + 1 < MATSLISE_ETA) {
            double dEta = -h * h * eta[i + 1] / 2;
            t.dt(0, 0) += dEta * hu[i];
            t.dt(0, 1) += dEta * hv[i];
            t.dt(1, 0) += dEta * hup[i];
            t.dt(1, 1) += dEta * hvp[i];
        }
    }
    delete[] eta;
    return t;
}

double Sector::prufer(double E, double delta, const Y &y0, const Y &y1) const {
    double theta0 = y0.theta();

    double theta1 = y1.theta();
    double ff = E - vs[0];
    if (ff > 0) {
        double f = sqrt(ff);
        double C = atan(y0.y[0] / y0.y[1] * f) / f;
        theta1 += round(((C + delta) * f - theta1) / M_PI) * M_PI;
    } else {
        if (y0.y[0] * y1.y[0] >= 0) {
            if (theta0 > 0 && theta1 < 0)
                theta1 += M_PI;
            else if (theta0 < 0 && theta1 > 0)
                theta1 -= M_PI;
        } else if (theta0 * theta1 > 0) {
            theta1 += M_PI;
        }
    }

    return theta1 - theta0;
}

Y Sector::propagate(double E, const Y &y0, double delta, double &theta) const {
    bool forward = delta >= 0;
    if (!forward)
        delta = -delta;
    const T &t = calculateT(E, delta);
    const Y y1 = forward ? t * y0 : t / y0;

    if (forward)
        theta += prufer(E, delta, y0, y1);
    else
        theta -= prufer(E, delta, y1, y0);

    return y1;
}

Sector::~Sector() {
    delete[]vs;
}
