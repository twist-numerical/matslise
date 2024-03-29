#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>
#include "../test.h"
#include "../../matslise/matslise.h"
#include "../../matslise/util/eigen.h"
#include "../../matslise/util/constants.h"
#include "test_problem.h"


using namespace matslise;
using namespace matslise::sector_builder;
using namespace std;
using namespace Eigen;

template<typename Scalar=double>
Scalar mathieu(Scalar x) {
    return 2 * cos(2 * x);
}

vector<double> correct = {-0.110248816992, 3.917024772998, 9.047739259809, 16.032970081406, 25.020840823290,
                          36.014289910628, 49.010418249424, 64.007937189250, 81.006250326633, 100.005050675159,
                          121.004166761269, 144.003496558993, 169.002976224509, 196.002564124230, 225.002232157099,
                          256.001960793939, 289.001736117776, 324.001547992331, 361.001388892288, 400.001253135325,
                          441.001136365493, 484.001035198089, 529.000946970769, 576.000869566047, 625.000801282700,
                          676.000740741253, 729.000686813595, 784.000638569932, 841.000595238360, 900.000556173742,
                          961.000520833511, 1024.000488758700, 1089.000459558945, 1156.000432900534, 1225.000408496817,
                          1296.000386100458, 1369.000365497137, 1444.000346500398, 1521.000328947413, 1600.000312695473,
                          1681.000297619080, 1764.000283607515, 1849.000270562795, 1936.000258397953, 2025.000247035591,
                          2116.000236406635, 2209.000226449289, 2304.000217108132, 2401.000208333344, 2500.000200080042,
                          2601.000192307701, 2704.000184979661, 2809.000178062684, 2916.000171526593, 3025.000165343920,
                          3136.000159489695, 3249.000153940890, 3364.000148676857, 3481.000143678164, 3600.000138927577,
                          3721.000134408604, 3844.000130106689, 3969.000126008066, 4096.000122100123, 4225.000118371213,
                          4356.000114810563, 4489.000111408200, 4624.000108154878, 4761.000105042017, 4900.000102061645,
                          5041.000099206349, 5184.000096469226, 5329.000093843843, 5476.000091324200, 5625.000088904693,
                          5776.000086580086, 5929.000084345478, 6084.000082196284, 6241.000080128204, 6400.000078137208,
                          6561.000076219510, 6724.000074371559, 6889.000072590010, 7056.000070871721, 7225.000069213730,
                          7396.000067613250, 7569.000066067651, 7744.000064574452, 7921.000063131311, 8100.000061736015,
                          8281.000060386470, 8464.000059080701, 8649.000057816833, 8836.000056593093, 9025.000055407798,
                          9216.000054259357, 9409.000053146274, 9604.000052067058, 9801.000051020470,
                          10000.000050004997, 10201.000049019618, 10404.000048063055, 10609.000047134234,
                          10816.000046232080, 11025.000045355584, 11236.000044503779, 11449.000043675756,
                          11664.000042870613, 11881.000042087557, 12100.000041325726, 12321.000040584424,
                          12544.000039862867, 12769.000039160412, 12996.000038476333, 13225.000037810047,
                          13456.000037160904, 13689.000036528345, 13924.000035911797, 14161.000035310730,
                          14400.000034724628, 14641.000034152999, 14884.000033595372, 15129.000033051290,
                          15376.000032520320, 15625.000032002043, 15876.000031496058, 16129.000031001979,
                          16384.000030519437, 16641.000030048071, 16900.000029587543, 17161.000029137522,
                          17424.000028697694, 17689.000028267747, 17956.000027847389, 18225.000027436341,
                          18496.000027034326, 18769.000026641086, 19044.000026256363, 19321.000025879912,
                          19600.000025511497, 19881.000025150901, 20164.000024797893, 20449.000024452263,
                          20736.000024113811, 21025.000023782337, 21316.000023457651, 21609.000023139572,
                          21904.000022827917, 22201.000022522516, 22500.000022223201, 22801.000021929820,
                          23104.000021642205, 23409.000021360211, 23716.000021083702, 24025.000020812513,
                          24336.000020546544, 24649.000020285614, 24964.000020029634, 25281.000019778472,
                          25600.000019533418, 25921.000019290117, 26244.000019054482, 26569.000018819625,
                          26896.000018590817, 27225.000018366140, 27556.000018145518, 27889.000017928851,
                          28224.000017716036, 28561.000017506991, 28900.000017301631, 29241.000017099854,
                          29584.000016901591, 29929.000016706752, 30276.000016515267, 30625.000016327052,
                          30976.000016142039, 31329.000015960151, 31684.000015781323, 32041.000015605485,
                          32400.000015432564, 32761.000015262507, 33124.000015095240, 33489.000014930709,
                          33856.000014768855, 34225.000014609621, 34596.000014452948, 34969.000014298777,
                          35344.000014147059, 35721.000013997749, 36100.000013850789, 36481.000013706129,
                          36864.000013563724, 37249.000013423531, 37636.000013285491, 38025.000013149576,
                          38416.000013015735, 38809.000012883931, 39204.000012754113, 39601.000012626253,
                          40000.000012500299, 40401.000012376229};

TEST_CASE("Solving the mathieu problem (first 200)", "[matslise][mathieu]") {
    Matslise<double> ms(&mathieu<>, 0, constants<double>::PI, 1e-10, UniformSectorBuilder<Matslise<>>(12));

    testProblem(ms, Y<>::Dirichlet(), Y<>::Dirichlet(), correct, 1e-7);
}

TEST_CASE("Solving the mathieu problem (first 200) (auto)", "[matslise][mathieu][auto]") {
    Matslise<double> ms(&mathieu<>, 0, constants<double>::PI, 1e-10);

    testProblem(ms, Y<>::Dirichlet(), Y<>::Dirichlet(), correct, 1e-7);
}

TEST_CASE("Solving the mathieu problem (first 200) (auto) (negative boundary conditions)",
          "[matslise][mathieu][auto]") {
    Matslise<double> ms(&mathieu<>, 0, constants<double>::PI, 1e-10);

    Y<double> y0({0, -1}, {0, 0});
    testProblem(ms, y0, y0, correct, 1e-7);
}

TEST_CASE("Mathieu problem eigenfunctions", "[mathieu][matslise][eigenfunctions]") {
    ArrayXd x(30);
    x << 0.000000000000000, 0.108330781158269, 0.216661562316537, 0.324992343474806, 0.433323124633075,
            0.541653905791344, 0.649984686949612, 0.758315468107881, 0.866646249266150, 0.974977030424419,
            1.083307811582687, 1.191638592740956, 1.299969373899225, 1.408300155057494, 1.516630936215762,
            1.624961717374031, 1.733292498532300, 1.841623279690568, 1.949954060848837, 2.058284842007106,
            2.166615623165375, 2.274946404323643, 2.383277185481912, 2.491607966640180, 2.599938747798450,
            2.708269528956718, 2.816600310114987, 2.924931091273256, 3.033261872431524, 3.141592653589793;

    vector<double> y0 = {-0.00000000000000, 0.05957648373348, 0.12058229947745, 0.18428570118601, 0.25163008589212,
                         0.32306954040196, 0.39841414637477, 0.47670540198060, 0.55615034750173, 0.63414569796321,
                         0.70741709118712, 0.77228194717948, 0.82501918024869, 0.86230069024582, 0.88161654418474,
                         0.88161654418474, 0.86230069024582, 0.82501918024869, 0.77228194717948, 0.70741709118712,
                         0.63414569796321, 0.55615034750173, 0.47670540198060, 0.39841414637477, 0.32306954040196,
                         0.25163008589212, 0.1842857016118601, 0.12058229947745, 0.05957648373348, 0.00000000000000};

    vector<double> dy0 = {0.54770139190616, 0.55442179455968, 0.57384238641313, 0.60372505966206, 0.64030804382969,
                          0.67833289106582, 0.71121342385303, 0.73143606944606, 0.73124805111186, 0.70362163948128,
                          0.64338934196303, 0.54835108763979, 0.42009484461026, 0.26427976960890, 0.09022152372952,
                          -0.09022152372952, -0.26427976960890, -0.42009484461026, -0.54835108763979, -0.64338934196303,
                          -0.70362163948128, -0.73124805111186, -0.73143606944606, -0.71121342385303, -0.67833289106582,
                          -0.64030804382969, -0.60372505966206, -0.57384238641313, -0.55442179455968,
                          -0.54770139190616};

    vector<double> y3 = {-0.00000000000000, 0.32417284029898, 0.59543237414050, 0.76858796800638, 0.81320927029333,
                         0.71916231964893, 0.49966578556769, 0.19086652240175, -0.15277306606069, -0.46727765184401,
                         -0.69110153718396, -0.77802995899298, -0.70806542774819, -0.49329149619403, -0.17665033540715,
                         0.17665033540715, 0.49329149619403, 0.70806542774819, 0.77802995899298, 0.69110153718396,
                         0.46727765184401, 0.15277306606069, -0.19086652240175, -0.49966578556769, -0.71916231964893,
                         -0.81320927029333, -0.76858796800638, -0.59543237414050, -0.32417284029898, -0.00000000000000};


    vector<double> dy3 = {3.07626133037988, 2.82599938325101, 2.11136070876314, 1.03774249508994, -0.23002968436917,
                          -1.48587759669341, -2.51014152163739, -3.10573501735402, -3.13821749458431, -2.57141674375545,
                          -1.48796775591086, -0.08518552873415, 1.35809781500210, 2.53962297948860, 3.20361729625403,
                          3.20361729625403, 2.53962297948859, 1.35809781500210, -0.08518552873415, -1.48796775591086,
                          -2.57141674375544, -3.13821749458430, -3.10573501735402, -2.51014152163740, -1.48587759669341,
                          -0.23002968436916, 1.03774249508994, 2.11136070876314, 2.82599938325101, 3.07626133037988};

    Matslise<double> ms(&mathieu<>, 0, constants<double>::PI, 1e-11, UniformSectorBuilder<Matslise<>>(8));
    Y<double> ystart({0, 1}, {0, 0});

    {
        auto eigenvalues = ms.eigenvaluesByIndex(0, 1, ystart, ystart);
        REQUIRE(eigenvalues.size() == 1);
        REQUIRE(0 == eigenvalues[0].first);
        double e = eigenvalues[0].second;

        REQUIRE_THAT(e, WithinAbs(-0.11024881635796, 1e-8));
        Array<double, Dynamic, 2> result = (*ms.eigenfunction(e, ystart, ystart))(x);
        REQUIRE(((unsigned long) result.rows()) == y0.size());
        for (Eigen::Index i = result.rows() - 1; i >= 0; --i)
            result.row(i) *= dy0[0] / result(0, 1);

        for (Eigen::Index i = 0; i < result.rows(); ++i) {
            REQUIRE_THAT(result(i, 0), WithinAbs(y0[i], 1e-8));
            REQUIRE_THAT(result(i, 1), WithinAbs(dy0[i], 1e-8));
        }
    }

    {
        auto eigenvalues = ms.eigenvaluesByIndex(3, 4, ystart, ystart);
        REQUIRE(eigenvalues.size() == 1);
        double e = eigenvalues.at(0).second;

        REQUIRE_THAT(e, WithinAbs(16.03297008140580, 1e-8));
        Array<double, Dynamic, 2> result = (*ms.eigenfunction(e, ystart, ystart))(x);
        REQUIRE(result.rows() == static_cast<long>(y0.size()));
        for (Eigen::Index i = result.rows() - 1; i >= 0; --i)
            result.row(i) *= dy3[0] / result(0, 1);

        for (Eigen::Index i = 0; i < result.rows(); ++i) {
            REQUIRE_THAT(result(i, 0), WithinAbs(y3[i], 1e-8));
            REQUIRE_THAT(result(i, 1), WithinAbs(dy3[i], 1e-8));
        }
    }

    {
        auto eigenpairs = ms.eigenpairsByIndex(3, 4, ystart, ystart);
        REQUIRE(eigenpairs.size() == 1);
        REQUIRE(get<0>(eigenpairs.at(0)) == 3);
        double e = get<1>(eigenpairs.at(0));

        REQUIRE_THAT(e, WithinAbs(16.03297008140580, 1e-8));
        Array<double, Dynamic, 2> result = (*get<2>(eigenpairs.at(0)))(x);
        REQUIRE(result.rows() == static_cast<long>(y0.size()));
        for (Eigen::Index i = result.rows() - 1; i >= 0; --i)
            result.row(i) *= dy3[0] / result(0, 1);

        for (Eigen::Index i = 0; i < result.rows(); ++i) {
            REQUIRE_THAT(result(i, 0), WithinAbs(y3[i], 1e-8));
            REQUIRE_THAT(result(i, 1), WithinAbs(dy3[i], 1e-8));
        }
    }
}

#ifdef MATSLISE_LONG_DOUBLE

TEST_CASE("Solving the mathieu problem (first 200) (auto) (long)", "[matslise][mathieu][auto][long]") {
    Matslise<long double> ms(&mathieu<long double>, 0, constants<long double>::PI, 1e-10);

    vector<long double> correct_ld(correct.size());
    transform(correct.begin(), correct.end(), correct_ld.begin(),
              [](double a) -> long double { return static_cast<long double>(a); });

    testProblem<long double>(ms, Y<long double>::Dirichlet(), Y<long double>::Dirichlet(), correct_ld, 1e-7);
}

#endif

#ifdef MATSLISE_QUADMATH

#include <boost/multiprecision/float128.hpp>

using boost::multiprecision::float128;

TEST_CASE("Solving the mathieu problem (first 200) (auto) (float128)", "[matslise][mathieu][auto][float128]") {
    Matslise<float128> ms(&mathieu<float128>, 0, constants<float128>::PI, 1e-10);

    vector<float128> correct_q(correct.size());
    transform(correct.begin(), correct.end(), correct_q.begin(),
              [](double a) -> float128 { return static_cast<float128>(a); });

    testProblem<float128>(ms, Y<float128>::Dirichlet(), Y<float128>::Dirichlet(), correct_q, 1e-7);
}

#endif
