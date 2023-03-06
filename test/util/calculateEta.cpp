#include <cmath>
#include <vector>
#include <tuple>
#include "../test.h"
#include "../../matslise/matslise.h"
#include "../../matslise/util/calculateEta.h"

using namespace std;
using namespace Eigen;

const int etaCount = 11;

template<typename Scalar>
void testEta(const Scalar &Z, const std::array<Scalar, etaCount> &correct, const Scalar &tolerance) {
    Array<Scalar, etaCount, 1> etas = calculateEta<Scalar, etaCount>(Z);

    for (int i = 0; i < etaCount; ++i)
        REQUIRE_THAT(etas[i], WithinAbs(correct[i], tolerance));
}

TEST_CASE("calculate eta", "[util][eta]") {
    testEta<double>(
            -50.0000000000000000, {
                    0.705347906308442312, 0.100248125275867068, -0.0121019956206515049, -0.00273108224275643165,
                    -0.0000310683118626130681, 0.0000502720811943628036, 9.67034085223756600e-6, 1.12203336360500845e-6,
                    9.83218574925508768e-8, 7.05588997566509405e-9, 4.32565441875114442e-10
            }, 1e-8);
    testEta<double>(
            -30.0000000000000000, {
                    0.692419111593747840, -0.131726455695091229, -0.0274715189096279690, 0.00164372996554024407,
                    0.00118967229124430631, 0.000222799202438996670, 0.0000271840176902221240, 2.54083307178155646e-6,
                    1.94893741431270334e-7, 1.27524349895832848e-8, 7.29921779721516930e-10
            }, 1e-8);
    testEta<double>(
            -10.0000000000000000, {
                    -0.999786072879325908, -0.00654070696893864021, 0.0993245365910387267, 0.0304514316742054820,
                    0.00529326217799886835, 0.000660140357178659639, 0.0000648001036609068404, 5.26607830913156049e-6,
                    3.65891435780344603e-7, 2.22293227573608555e-8, 1.20070510947899410e-9
            }, 1e-8);
    testEta<double>(
            -1.00000000000000000, {
                    0.540302305868139717, 0.841470984807896507, 0.301168678939756789, 0.0620350520113738611,
                    0.00900658111711251626, 0.00101101580841375271, 0.0000925611586112581636, 7.15693631008708557e-6,
                    4.79013419873948858e-7, 2.82649880221472943e-8, 1.49137650255514566e-9
            }, 1e-8);
    testEta<double>(
            -0.200000000000000000, {
                    0.901655595150045281, 0.966998417099577368, 0.326714109747660432, 0.0657195607170196462,
                    0.00941846918718899323, 0.00104861796651653215, 0.0000954625572989805510, 7.35081886126956194e-6,
                    4.90439487618770872e-7, 2.88672650600056869e-8, 1.52009200662902872e-9
            }, 1e-8);
    testEta<double>(
            -0.0100000000000000000, {
                    0.995004165278025766, 0.998334166468281523, 0.333000119025575697, 0.0666190608445568706,
                    0.00951851972086556705, 0.00105772015020987319, 0.0000961631023291644604, 7.39754109358770711e-6,
                    4.93188747573197342e-7, 2.90120010253018994e-8, 1.52698569349481951e-9
            }, 1e-8);
    testEta<double>(
            0.0100000000000000000, {
                    1.00500416805580360, 1.00166750019844026, 0.333666785736334075, 0.0667142989438033032,
                    0.00952910173175591127, 0.00105868215119242973, 0.0000962371024043736860, 7.40247443191801826e-6,
                    4.93478943944855538e-7, 2.90272745185193289e-8, 1.52771300269464561e-9
            }, 1e-8);
    testEta<double>(
            0.200000000000000000, {
                    1.10167781754863464, 1.03366825838545201, 0.340047795815913168, 0.0676243546885625067,
                    0.00963011186550317482, 0.00106785815020141484, 0.0000969425684522064404, 7.44948613571995716e-6,
                    4.96243439234986595e-7, 2.91727359757911759e-8, 1.53463823268302235e-9
            }, 1e-8);
    testEta<double>(
            1.00000000000000000, {
                    1.54308063481524378, 1.17520119364380146, 0.367879441171442322, 0.0715628701294744921,
                    0.0100650905240698611, 0.00110723646098546428, 0.0000999623752006825919, 7.65033377795577012e-6,
                    5.08036087257580248e-7, 2.97924690920664063e-8, 1.56411269245133990e-9
            }, 1e-8);
    testEta<double>(
            2.00000000000000000, {
                    2.17818355660857086, 1.36829887200859068, 0.404942342299990092, 0.0767359225543102008,
                    0.0106313647642195443, 0.00115818460238669522, 0.000103851671369643672, 7.90810866030741339e-6,
                    5.23129392823649059e-7, 3.05838839763387562e-8, 1.60168261294510167e-9
            }, 1e-8);
}


#ifdef MATSLISE_LONG_DOUBLE

TEST_CASE("calculate eta (long)", "[util][eta][long]") {
    testEta<long double>(
            -50.000000000000000000000l, {
                    0.70534790630844231151456l, 0.10024812527586706813744l, -0.012101995620651504867543l,
                    -0.0027310822427564316548013l, -0.000031068311862613068129280l, 0.000050272081194362803557927l,
                    9.6703408522375660030124e-6l, 1.1220333636050084495042e-6l, 9.8321857492550876810843e-8l,
                    7.0558899756650940531688e-9l, 4.3256544187511444186056e-10l
            }, 1e-8l);
    testEta<long double>(
            -30.000000000000000000000l, {
                    0.69241911159374784000566l, -0.13172645569509122915402l, -0.027471518909627968971990l,
                    0.0016437299655402440746018l, 0.0011896722912443063115000l, 0.00022279920243899667019660l,
                    0.000027184017690222124008980l, 2.5408330717815564634062e-6l, 1.9489374143127033384336e-7l,
                    1.2752434989583284808140e-8l, 7.2992177972151692983385e-10l
            }, 1e-8l);
    testEta<long double>(
            -10.000000000000000000000l, {
                    -0.99978607287932590757522l, -0.0065407069689386402127975l, 0.099324536591038726736242l,
                    0.030451431674205482042152l, 0.0052932621779988683474520l, 0.00066014035717865963900115l,
                    0.000064800103660906840355835l, 5.2660783091315604913035e-6l, 3.6589143578034460311107e-7l,
                    2.2229322757360855536252e-8l, 1.2007051094789941005210e-9l
            }, 1e-8l);
    testEta<long double>(
            -1.0000000000000000000000l, {
                    0.54030230586813971740094l, 0.84147098480789650665250l, 0.30116867893975678925157l,
                    0.062035052011373861102195l, 0.0090065811171125162594084l, 0.0010110158084137527136639l,
                    0.000092561158611258163566821l, 7.1569363100870855711166e-6l, 4.7901341987394885769542e-7l,
                    2.8264988022147294314731e-8l, 1.4913765025551456549984e-9l
            }, 1e-8l);
    testEta<long double>(
            -0.20000000000000000000000l, {
                    0.90165559515004528111994l, 0.96699841709957736757398l, 0.32671410974766043227021l,
                    0.065719560717019646183259l, 0.0094184691871889932304359l, 0.0010486179665165321489603l,
                    0.000095462557298980551031640l, 7.3508188612695619389223e-6l, 4.9043948761877087175364e-7l,
                    2.8867265060005686911638e-8l, 1.5200920066290287209913e-9l
            }, 1e-8l);
    testEta<long double>(
            -0.010000000000000000000000l, {
                    0.99500416527802576609556l, 0.99833416646828152306814l, 0.33300011902557569725800l,
                    0.066619060844556870585691l, 0.0095185197208655670453669l, 0.0010577201502098731877746l,
                    0.000096163102329164460440516l, 7.3975410935877071087642e-6l, 4.9318874757319734185056e-7l,
                    2.9012001025301899414453e-8l, 1.5269856934948195148215e-9l
            }, 1e-8l);
    testEta<long double>(
            0.010000000000000000000000l, {
                    1.0050041680558035989880l, 1.0016675001984402582373l, 0.33366678573633407506846l,
                    0.066714298943803303191151l, 0.0095291017317559112705559l, 0.0010586821511924297259662l,
                    0.000096237102404373685980327l, 7.4024744319180182634575e-6l, 4.9347894394485553801122e-7l,
                    2.9027274518519328915005e-8l, 1.5277130026946456136098e-9l
            }, 1e-8l);
    testEta<long double>(
            0.20000000000000000000000l, {
                    1.1016778175486346400269l, 1.0336682583854520063587l, 0.34004779581591316834119l,
                    0.067624354688562506675618l, 0.0096301118655031748154775l, 0.0010678581502014148363764l,
                    0.000096942568452206440449421l, 7.4494861357199571638747e-6l, 4.9624343923498659524594e-7l,
                    2.9172735975791175928006e-8l, 1.5346382326830223491700e-9l
            }, 1e-8l);
    testEta<long double>(
            1.0000000000000000000000l, {
                    1.5430806348152437784779l, 1.1752011936438014568824l, 0.36787944117144232159552l,
                    0.071562870129474492095811l, 0.010065090524069861116471l, 0.0011072364609854642805131l,
                    0.000099962375200682591853594l, 7.6503337779557701235217e-6l, 5.0803608725758024781159e-7l,
                    2.9792469092066406347822e-8l, 1.5641126924513398986164e-9l
            }, 1e-8l);
    testEta<long double>(
            2.0000000000000000000000l, {
                    2.1781835566085708639892l, 1.3682988720085906790060l, 0.40494234229999009249163l,
                    0.076735922554310200765535l, 0.010631364764219544331978l, 0.0011581846023866952208447l,
                    0.00010385167136964367218769l, 7.9081086603074133900515e-6l, 5.2312939282364905851193e-7l,
                    3.0583883976338756186291e-8l, 1.6016826129451016724897e-9l
            }, 1e-8l);
}

#endif

#ifdef MATSLISE_QUADMATH

#include <boost/multiprecision/float128.hpp>

using boost::multiprecision::float128;

TEST_CASE("calculate eta (float128)", "[util][eta][float128]") {
    testEta<float128>(
            -50.000000000000000000000000000000000000q, {
                    0.70534790630844231151456449684134073842q, 0.10024812527586706813743756939417117421q,
                    -0.012101995620651504867542538548943391284q, -0.0027310822427564316548013037008200269612q,
                    -0.000031068311862613068129279599103134870434q, 0.000050272081194362803557926930141961657363q,
                    9.6703408522375660030124394076157957340e-6q, 1.1220333636050084495041980668362419142e-6q,
                    9.8321857492550876810842709225106983018e-8q, 7.0558899756650940531688514308072566209e-9q,
                    4.3256544187511444186055530197232759075e-10q
            }, 1e-24q);
    testEta<float128>(
            -30.000000000000000000000000000000000000q, {
                    0.69241911159374784000566456437372974808q, -0.13172645569509122915402275772007756507q,
                    -0.027471518909627968971989577403126910438q, 0.0016437299655402440746018008503565611251q,
                    0.0011896722912443063114999527218303238688q, 0.00022279920243899667019659560674852353188q,
                    0.000027184017690222124008980257963546263938q, 2.5408330717815564634062410283495123812e-6q,
                    1.9489374143127033384336251349991323394e-7q, 1.2752434989583284808139889138306204262e-8q,
                    7.2992177972151692983385339504307461695e-10q
            }, 1e-24q);
    testEta<float128>(
            -10.000000000000000000000000000000000000q, {
                    -0.99978607287932590757522142516383587566q, -0.0065407069689386402127974939033799063819q,
                    0.099324536591038726736242393126045596928q, 0.030451431674205482042152467328151669717q,
                    0.0052932621779988683474519943514712751655q, 0.00066014035717865963900114931321472564419q,
                    0.000064800103660906840355834946746125563219q, 5.2660783091315604913035100992655551228e-6q,
                    3.6589143578034460311106845443266533762e-7q, 2.2229322757360855536251671722442494161e-8q,
                    1.2007051094789941005209964848857063112e-9q
            }, 1e-24q);
    testEta<float128>(
            -1.0000000000000000000000000000000000000q, {
                    0.54030230586813971740093660744297660373q, 0.84147098480789650665250232163029899962q,
                    0.30116867893975678925156571418732239589q, 0.062035052011373861102194820931668188048q,
                    0.0090065811171125162594083904710185443507q, 0.0010110158084137527136639123654616224069q,
                    0.000092561158611258163566820818136057310989q, 7.1569363100870855711166340350080140233e-6q,
                    4.7901341987394885769542431904687131323e-7q, 2.8264988022147294314730750695055675278e-8q,
                    1.4913765025551456549984427690751664944e-9q
            }, 1e-24q);
    testEta<float128>(
            -0.20000000000000000000000000000000000000q, {
                    0.90165559515004528111993664148377875958q, 0.96699841709957736757397867981290111066q,
                    0.32671410974766043227021019164561175542q, 0.065719560717019646183259475619670777926q,
                    0.0094184691871889932304359322637106710688q, 0.0010486179665165321489602511315195977783q,
                    0.000095462557298980551031639599828544678900q, 7.3508188612695619389223329719684480776e-6q,
                    4.9043948761877087175364403522573053978e-7q, 2.8867265060005686911637782087550096015e-8q,
                    1.5200920066290287209913013131054623206e-9q
            }, 1e-24q);
    testEta<float128>(
            -0.010000000000000000000000000000000000000q, {
                    0.99500416527802576609556198780387029484q, 0.99833416646828152306814198410622026990q,
                    0.33300011902557569725799963023499750606q, 0.066619060844556870585690659877224827414q,
                    0.0095185197208655670453669151126631013509q, 0.0010577201502098731877745911416882042146q,
                    0.000096163102329164460440516253073658062545q, 7.3975410935877071087642122034473378755e-6q,
                    4.9318874757319734185055711573298359054e-7q, 2.9012001025301899414453254741598258438e-8q,
                    1.5269856934948195148214874186802905265e-9q
            }, 1e-24q);
    testEta<float128>(
            0.010000000000000000000000000000000000000q, {
                    1.0050041680558035989879784429683416447q, 1.0016675001984402582372938352190502352q,
                    0.33366678573633407506846077492914095604q, 0.066714298943803303191151043162736702408q,
                    0.0095291017317559112705559115457444003487q, 0.0010586821511924297259662342525899966381q,
                    0.000096237102404373685980327243443060626389q, 7.4024744319180182634574716329747776974e-6q,
                    4.9347894394485553801122143885163225569e-7q, 2.9027274518519328915005020029386202350e-8q,
                    1.5277130026946456136098352066815745954e-9q
            }, 1e-24q);
    testEta<float128>(
            0.20000000000000000000000000000000000000q, {
                    1.1016778175486346400269158064791200840q, 1.0336682583854520063586787924895828541q,
                    0.34004779581591316834118506994768614930q, 0.067624354688562506675617913232622031218q,
                    0.0096301118655031748154775189228799660428q, 0.0010678581502014148363764038623113445939q,
                    0.000096942568452206440449420810389323487987q, 7.4494861357199571638747401439311302529e-6q,
                    4.9624343923498659524594259109397349832e-7q, 2.9172735975791175928006387607638890719e-8q,
                    1.5346382326830223491700088205617805024e-9q
            }, 1e-24q);
    testEta<float128>(
            1.0000000000000000000000000000000000000q, {
                    1.5430806348152437784779056207570616826q, 1.1752011936438014568823818505956008152q,
                    0.36787944117144232159552377016146086745q, 0.071562870129474492095810540111218212818q,
                    0.010065090524069861116471069605369803354q, 0.0011072364609854642805130528736295893376q,
                    0.000099962375200682591853593742703499316283q, 7.6503337779557701235217038910968584558e-6q,
                    5.0803608725758024781159211924015635682e-7q, 2.9792469092066406347822102494513103611e-8q,
                    1.5641126924513398986163768334335954292e-9q
            }, 1e-24q);
    testEta<float128>(
            2.0000000000000000000000000000000000000q, {
                    2.1781835566085708639892220678201252834q, 1.3682988720085906790059611763827516729q,
                    0.40494234229999009249163044571868680526q, 0.076735922554310200765534919613345628577q,
                    0.010631364764219544331977923825979331186q, 0.0011581846023866952208447264157451551375q,
                    0.00010385167136964367218769304213646747425q, 7.9081086603074133900514761220064603558e-6q,
                    5.2312939282364905851192627519174481369e-7q, 3.0583883976338756186290997065144075206e-8q,
                    1.6016826129451016724896625421477675952e-9q
            }, 1e-24q);
}

#endif
