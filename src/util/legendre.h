//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_LEGENDRE_H
#define SCHRODINGER_LEGENDRE_H

#include <functional>
#include <iostream>

namespace legendre {
    template<class D = double>
    class array2d {
        int rows, cols;
        D *values;
    public:
        array2d(int rows, int cols, D *values) : rows(rows), cols(cols), values(values) {}

        D *operator[](int row) const {
            return values + cols * row;
        }

        virtual ~array2d() {
            delete[] values;
        }
    };

    struct LegendreData {
        int n;
        array2d<double> polynomial;
        double *nodes;
        double *weights;
    };

    static LegendreData legendreData[] = {
            {
                    0,
                    array2d<double>(0,0,new double[0]),
                    new double[0],
                    new double[0]
            },
            {
                    2,
                    array2d<double>(2, 2, new double[4]{
                            1.0,1.0,
                            -0.5773502691896257,0.5773502691896256
                    }),
                    new double[2]{-0.5773502691896258,0.5773502691896256},
                    new double[2]{0.9999999999999998,1.0000000000000004}
            },
            {
                    4,
                    array2d<double>(4, 4, new double[16]{
                            1.0,1.0,1.0,1.0,
                            -0.8611363115940536,-0.3399810435848565,0.33998104358485626,0.8611363115940531,
                            0.6123336207187163,-0.32661933500442786,-0.3266193350044281,0.6123336207187152,
                            -0.30474698495521024,0.4117279996728997,-0.4117279996728996,0.3047469849552084
                    }),
                    new double[4]{-0.8611363115940536,-0.33998104358485653,0.3399810435848563,0.8611363115940531},
                    new double[4]{0.3478548451374517,0.6521451548625459,0.6521451548625463,0.3478548451374526}
            },
            {
                    6,
                    array2d<double>(6, 6, new double[36]{
                            1.0,1.0,1.0,1.0,1.0,1.0,
                            -0.9324695142031534,-0.6612093864662638,-0.23861918608319677,0.23861918608319688,0.6612093864662635,0.9324695142031517,
                            0.8042490923773973,0.15579677912663953,-0.41459132604948906,-0.414591326049489,0.15579677912663886,0.8042490923773926,
                            -0.6282499246436956,0.26911576974460033,0.32396186535393506,-0.3239618653539352,-0.2691157697446009,0.6282499246436872,
                            0.42200500927063245,-0.4282458620971208,0.175662340429804,0.17566234042980383,-0.4282458620971208,0.42200500927062046,
                            -0.2057123110596345,0.2943957149254359,-0.3346190207410407,0.3346190207410408,-0.29439571492543526,0.20571231105961982
                    }),
                    new double[6]{-0.9324695142031534,-0.6612093864662637,-0.23861918608319677,0.23861918608319688,0.6612093864662635,0.9324695142031517},
                    new double[6]{0.1713244923791671,0.3607615730481392,0.46791393457269126,0.46791393457269126,0.3607615730481398,0.1713244923791712}
            },
            {
                    8,
                    array2d<double>(8, 8, new double[64]{
                            1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                            -0.9602898564975371,-0.7966664774136296,-0.5255324099163281,-0.1834346424956499,0.18343464249564978,0.5255324099163294,0.7966664774136274,0.9602898564975355,
                            0.8832349127380905,0.4520162143519618,-0.08572352919130473,-0.4495275978987397,-0.44952759789873975,-0.08572352919130263,0.45201621435195655,0.883234912738086,
                            -0.7734093083464346,-0.06906595709361621,0.42543909474828384,0.25972131868457254,-0.2597213186845724,-0.4254390947482831,0.06906595709360897,0.7734093083464262,
                            0.6372937644666821,-0.2427227284567702,-0.3269759103939731,0.25377239575159893,0.2537723957515991,-0.326975910393975,-0.2427227284567766,0.6372937644666692,
                            -0.4828486810505258,0.4033170755970771,-0.031045687085553008,-0.2915682225895845,0.29156822258958437,0.031045687085549847,-0.40331707559707947,0.4828486810505085,
                            0.31899212911049485,-0.3867979150966231,0.3023917023728722,-0.11342352322434288,-0.11342352322434308,0.30239170237287083,-0.38679791509661965,0.3189921291104737,
                            -0.1550188128903624,0.2265762384000079,-0.26852031408771393,0.2885549685956877,-0.28855496859568763,0.268520314087716,-0.22657623839999913,0.15501881289033861
                    }),
                    new double[8]{-0.9602898564975371,-0.7966664774136296,-0.525532409916328,-0.18343464249564984,0.18343464249564978,0.5255324099163294,0.7966664774136274,0.9602898564975355},
                    new double[8]{0.10122853629037425,0.22238103445337207,0.31370664587788794,0.36268378337836193,0.36268378337836193,0.3137066458778873,0.22238103445337357,0.10122853629037788}
            },
            {
                    10,
                    array2d<double>(10, 10, new double[100]{
                            1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                            -0.9739065285171686,-0.8650633666889882,-0.6794095682990249,-0.433395394129247,-0.14887433898163116,0.14887433898163116,0.43339539412924677,0.6794095682990237,0.8650633666889889,0.9739065285171682,
                            0.922740889432544,0.6225019425809303,0.19239604224440102,-0.218252648521332,-0.46675464678917356,-0.46675464678917356,-0.21825264852133228,0.19239604224439855,0.622501942580932,0.9227408894325426,
                            -0.8485012749020432,-0.32079713257316583,0.2350801921931675,0.44657975046225584,0.21506254183332568,-0.21506254183332568,-0.4465797504622558,-0.2350801921931699,0.32079713257316855,0.8485012749020407,
                            0.7540759623195132,0.018765776238156177,-0.4237995624971213,-0.17501532579202875,0.2940357210203751,0.2940357210203751,-0.17501532579202833,-0.42379956249712175,0.0187657762381594,0.7540759623195094,
                            -0.6431180849398782,0.22741725203053179,0.33021610628813863,-0.2207322953892739,-0.25084390595367273,0.25084390595367273,0.22073229538927427,-0.33021610628813636,-0.22741725203052893,0.6431180849398727,
                            0.5198876842061506,-0.37631042528706177,-0.05814566531984909,0.32123076511505144,-0.17656536292520286,-0.17656536292520286,0.32123076511505133,-0.058145665319845136,-0.3763104252870602,0.5198876842061435,
                            -0.3890682310047542,0.4096310303233839,-0.20967646569634207,-0.06935219576565162,0.26382601538929634,-0.26382601538929634,0.06935219576565106,0.20967646569634527,-0.40963103032338427,0.3890682310047457,
                            0.2555659454711607,-0.3351473744834937,0.31798282660714994,-0.22472019031770138,0.08085046072097912,0.08085046072097912,-0.22472019031770174,0.31798282660715005,-0.33514737448349624,0.255565945471151,
                            -0.12430099765549005,0.18351721458258377,-0.22169756095640109,0.24561037653348733,-0.25724773603885603,0.25724773603885603,-0.24561037653348702,0.22169756095639767,-0.18351721458258796,0.12430099765547961
                    }),
                    new double[10]{-0.9739065285171685,-0.8650633666889883,-0.679409568299025,-0.433395394129247,-0.14887433898163116,0.1488743389816312,0.4333953941292468,0.6794095682990237,0.8650633666889889,0.9739065285171682},
                    new double[10]{0.06667134430869694,0.14945134915057698,0.21908636251598265,0.26926671930999646,0.2955242247147528,0.295524224714753,0.2692667193099972,0.21908636251598257,0.1494513491505774,0.06667134430869787}
            },
            {
                    12,
                    array2d<double>(12, 12, new double[144]{
                            1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                            -0.9815606342467325,-0.9041172563704614,-0.7699026741943014,-0.5873179542866223,-0.367831498998179,-0.12523340851146902,0.1252334085114689,0.3678314989981788,0.5873179542866254,0.7699026741942895,0.9041172563704775,0.9815606342467244,
                            0.9451919180542715,0.7261420199002758,0.38912519159730485,0.017413569141134365,-0.29704998251712894,-0.47647489008889926,-0.4764748900888993,-0.29704998251712916,0.017413569141139843,0.3891251915972774,0.7261420199003195,0.9451919180542476,
                            -0.8918982081192389,-0.491451047032263,0.013954240117642786,0.3744997998497487,0.42732823324321667,0.1829398966009134,-0.18293989660091323,-0.42732823324321656,-0.3744997998497453,-0.013954240117677779,0.49145104703233755,0.8918982081191925,
                            0.823147360438124,0.23296988657039347,-0.3106448555680931,-0.3979734754307326,-0.05228588615882107,0.3172633406595553,0.3172633406595554,-0.05228588615882065,-0.39797347543073525,-0.31064485556811944,0.2329698865704925,0.8231473604380495,
                            -0.7408257148469124,0.014023067232479284,0.4193359569537251,0.12112670153076838,-0.3072442740065282,-0.21786946246448732,0.21786946246448716,0.30724427400652843,-0.12112670153077612,-0.4193359569537269,-0.014023067232370985,0.7408257148468059,
                            0.6471803569425079,-0.2173854834408953,-0.3330170572081713,0.2012214375126257,0.25076412855413655,-0.214364468992149,-0.21436446899214917,0.2507641285541362,0.20122143751261887,-0.3330170572081428,-0.21738548344079875,0.6471803569423674,
                            -0.5447505160201739,0.3529867379755224,0.11672337941959854,-0.32330181838364047,0.09205133644652115,0.23660135504145957,-0.23660135504145943,-0.09205133644652165,0.32330181838364086,-0.11672337941954888,-0.35298673797545954,0.5447505160199992,
                            0.4362903039459244,-0.4081778290021945,0.12289184638937803,0.17995804703312407,-0.28290495199094995,0.13201192133636622,0.13201192133636644,-0.28290495199094995,0.1799580470331323,0.12289184638942735,-0.4081778290021832,0.43629030394571766,
                            -0.3246852731783611,0.3833107352182204,-0.2824708861420405,0.08773783137485011,0.11473692244778817,-0.24153999879715496,0.24153999879715488,-0.11473692244778762,-0.08773783137484023,0.2824708861420654,-0.3833107352182693,0.32468527317812634,
                            0.21286546352456684,-0.2910998693974244,0.3026000104328627,-0.2598692492393927,0.17442713386048023,-0.06133786225440428,-0.06133786225440457,0.1744271338604807,-0.2598692492393896,0.3026000104328483,-0.2910998693975303,0.21286546352431013,
                            -0.1037177104846144,0.15398630618683398,-0.1879740764266154,0.21161500718203535,-0.22679317280626068,0.23424659352281785,-0.2342465935228178,0.22679317280626043,-0.21161500718204243,0.18797407642656475,-0.15398630618698128,0.10371771048434342
                    }),
                    new double[12]{-0.9815606342467326,-0.9041172563704614,-0.7699026741943014,-0.5873179542866224,-0.36783149899817896,-0.125233408511469,0.12523340851146894,0.3678314989981788,0.5873179542866254,0.7699026741942895,0.9041172563704775,0.9815606342467244},
                    new double[12]{0.04717533638647901,0.10693932599533276,0.16007832854333748,0.20316742672306162,0.23349253653835503,0.24914704581340288,0.24914704581340272,0.23349253653835528,0.20316742672306173,0.1600783285433453,0.10693932599531966,0.04717533638649818}
            },
            {
                    14,
                    array2d<double>(14, 14, new double[196]{
                            1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                            -0.9862838086967527,-0.9284348836637202,-0.8272013150696356,-0.6872929048117413,-0.5152486363581417,-0.3191123689278914,-0.10805494870734345,0.10805494870734345,0.3191123689278912,0.515248636358144,0.6872929048117465,0.8272013150696276,0.9284348836637257,0.986283808696745,
                            0.959133626946059,0.7929869998054985,0.5263930234794018,0.20855730550684198,-0.10177826409661324,-0.347250943995844,-0.48248619208978005,-0.48248619208978005,-0.3472509439958442,-0.10177826409660964,0.20855730550685272,0.5263930234793819,0.7929869998055139,0.9591336269460363,
                            -0.9191074052579806,-0.6081047324096206,-0.17425412539631407,0.21929534267198222,0.43090094388322575,0.39742836487019484,0.15892833370200005,-0.15892833370200005,-0.3974283648701947,-0.4309009438832246,-0.21929534267197157,0.17425412539628504,0.6081047324096481,0.9191074052579362,
                            0.8670260962441098,0.39328463150352283,-0.142544094662306,-0.420178212014379,-0.3122032684747594,0.03849567075936706,0.331811906411287,0.331811906411287,0.03849567075936747,-0.3122032684747628,-0.4201782120143763,-0.14254409466233559,0.3932846315035619,0.8670260962440378,
                            -0.8039549165912276,-0.17076672204636759,0.3516460929251775,0.3443776328355535,-0.055168600079866355,-0.34005469233723745,-0.19167972031156347,0.19167972031156347,0.34005469233723756,0.05516860007986095,-0.34437763283556255,-0.3516460929251963,0.1707667220464148,0.8039549165911235,
                            0.7311790679160481,-0.03706992643888643,-0.4144971237105562,-0.08378004663153875,0.31228289132402065,0.1668656481679595,-0.23853802770748195,-0.23853802770748195,0.1668656481679591,0.31228289132401865,-0.08378004663155576,-0.41449712371055486,-0.03706992643883691,0.7311790679159085,
                            -0.6501744982455993,0.21028878560106257,0.33535239974090436,-0.18824388371754666,-0.25153310576937626,0.19258465062701025,0.21216515834317964,-0.21216515834317964,-0.19258465062701063,0.2515331057693803,0.1882438837175319,-0.33535239974088005,-0.21028878560101805,0.6501744982454222,
                            0.5625744039119012,-0.333637772232176,-0.15744741564845247,0.31589257640248947,-0.0302436116336443,-0.26123771229834375,0.16573547055231358,0.16573547055231358,-0.2612377122983435,-0.030243611633637518,0.3158925764024872,-0.15744741564841092,-0.33363777223214397,0.5625744039116867,
                            -0.47013227244992273,0.3981806446012169,-0.05208079335434146,-0.24277014220109658,0.25301950002803086,-0.013720672910446175,-0.22241848986969306,0.22241848986969306,0.013720672910446766,-0.25301950002802803,0.24277014220110985,0.05208079335438715,-0.3981806446012038,0.4701322724496725,
                            0.37468234817931045,-0.40212712584124577,0.2235571455135115,0.03271965408416672,-0.22047985901649736,0.24343297029641012,-0.10349842831965672,-0.10349842831965672,0.24343297029641028,-0.22047985901650177,0.03271965408418851,0.2235571455135451,-0.40212712584125576,0.3746823481790276,
                            -0.27809846156265167,0.35077449375046715,-0.3056958296780352,0.17776851944656158,-0.01314128359308286,-0.13582956175007163,0.22354897848433866,-0.22354897848433866,0.13582956175007108,0.013141283593074926,-0.17776851944654473,0.3056958296780434,-0.35077449375050096,0.27809846156234114,
                            0.1822521997391779,-0.25558674754890404,0.27974393522781477,-0.2641695136299633,0.21508434196350806,-0.1400691774324125,0.04857537693591887,0.04857537693591887,-0.14006917743241304,0.21508434196350434,-0.2641695136299629,0.27974393522779223,-0.2555867475489588,0.18225219973884607,
                            -0.08897140798470174,0.13254595235832906,-0.16282798626489697,0.1850633520213324,-0.20098864956539614,0.2113384551236013,-0.21644676833794835,0.21644676833794835,-0.21133845512360105,0.20098864956540072,-0.18506335202135,0.16282798626484923,-0.13254595235839836,0.08897140798435631
                    }),
                    new double[14]{-0.9862838086967527,-0.9284348836637202,-0.8272013150696356,-0.6872929048117412,-0.5152486363581417,-0.3191123689278915,-0.10805494870734346,0.10805494870734345,0.3191123689278912,0.515248636358144,0.6872929048117465,0.8272013150696276,0.9284348836637257,0.986283808696745},
                    new double[14]{0.035119460331904884,0.08015808715957283,0.1215185706879927,0.15720316715815244,0.1855383974779421,0.20519846372129524,0.21526385346315782,0.21526385346315782,0.2051984637212956,0.1855383974779424,0.15720316715816968,0.12151857068796204,0.0801580871596417,0.03511946033192289}
            },
            {
                    16,
                    array2d<double>(16, 16, new double[256]{
                            1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                            -0.9894009349916444,-0.9445750230731487,-0.865631202388051,-0.7554044083547979,-0.6178762444027379,-0.45801677765720394,-0.2816035507792618,-0.09501250983763732,0.09501250983763732,0.2816035507792616,0.4580167776571993,0.6178762444027687,0.7554044083547129,0.8656312023881924,0.944575023072992,0.9894009349917203,
                            0.9683713152435102,0.838332961320459,0.6239760678216745,0.3559537302427935,0.07265658009584788,-0.18533094707676712,-0.3810491602827676,-0.4864589344615293,-0.4864589344615293,-0.38104916028276775,-0.18533094707677353,0.07265658009590488,0.35595373024260074,0.6239760678220414,0.8383329613200151,0.9683713152437355,
                            -0.9372451845405994,-0.6900639450883236,-0.3231344548243728,0.055454577243630704,0.3370962048672099,0.44681899040523376,0.3665770281146469,0.14037448038531258,-0.14037448038531258,-0.3665770281146467,-0.4468189904052341,-0.337096204867168,-0.0554545772438671,0.3231344548249551,0.6900639450875103,0.937245184541043,
                            0.8965162218939167,0.5119303210037204,0.0215196657933453,-0.3402739038803531,-0.41898897493697923,-0.21914032950996268,0.10513543289748646,0.3415038703736622,0.3415038703736622,0.10513543289748684,-0.2191403295099545,-0.4189889749369948,-0.3402739038805128,0.021519665794032004,0.5119303210025198,0.8965162218946403,
                            -0.8468290310861524,-0.31835071452242164,0.22497699434326976,0.41831627087586337,0.19631303781061657,-0.17678928668578323,-0.3465533426816821,-0.1707044360264939,0.1707044360264939,0.3465533426816821,0.17678928668579208,-0.1963130378106906,-0.41831627087583934,-0.2249769943426601,0.3183507145208866,0.8468290310872088,
                            0.788967779502033,0.1246859772750776,-0.3749694160533141,-0.2957679977943452,0.1267793478096692,0.33106645018055536,0.09130333428033192,-0.2548516209813433,-0.2548516209813433,0.09130333428033145,0.3310664501805545,0.12677934780958725,-0.2957679977941136,-0.3749694160529771,0.12468597727332842,0.7889677795034641,
                            -0.7238423966873992,0.05414655841928159,0.4099637111630092,0.05637431667691234,-0.31374564880939554,-0.13007230471609918,0.24929608505469097,0.19128711629316048,-0.19128711629316048,-0.24929608505469125,0.13007230471608797,0.3137456488093722,-0.05637431667656131,-0.4099637111630885,-0.05414655842107056,0.7238423966892346,
                            0.6524725880652941,-0.20499801761861386,-0.3372968488836689,0.17894935931556022,0.2525480391366326,-0.17797946040452278,-0.2115204101454839,0.18891766394876813,0.18891766394876813,-0.2115204101454835,-0.17797946040453277,0.2525480391366953,0.17894935931586373,-0.3372968488842011,-0.2049980176202358,0.6524725880675496,
                            -0.5759699593556162,0.31762662838007,0.1870955352589705,-0.30544886963359935,-0.015863687414016118,0.2695974755300889,-0.10908504499224544,-0.20393768158007694,0.20393768158007694,0.10908504499224596,-0.269597475530086,0.015863687414124795,0.30544886963369167,-0.1870955352598602,-0.317626628381313,0.5759699593582938,
                            0.49551858173697294,-0.38554392582167446,-0.004148728985249368,0.2773466796475966,-0.20866982357964192,-0.07443080294886509,0.248733967543726,-0.1332102986968465,-0.1332102986968465,0.2487339675437259,-0.07443080294885118,-0.2086698235795699,0.27734667964740656,-0.00414872898628394,-0.38554392582235103,0.4955185817400595,
                            -0.412354356004184,0.40649201194567547,-0.16323079073928637,-0.12229075427996129,0.2605646853931957,-0.1800066425703874,-0.034552844339255505,0.20956057790221388,-0.20956057790221388,0.03455284433925487,0.18000664257039758,-0.2605646853932218,0.1222907542795583,0.16323079073838437,-0.40649201194565016,0.4123543560076518,
                            0.32774355538307415,-0.3825122876585194,0.27462352740741347,-0.0771747544824076,-0.1172963927518916,0.2262498556115064,-0.2093565799087887,0.08394692723887211,0.08394692723887211,-0.20935657990878898,0.226249855611501,-0.11729639275200388,-0.07717475448283682,0.2746235274069095,-0.3825122876577314,0.32774355538688094,
                            -0.24296094086825196,0.3196065138593935,-0.3064844513055187,0.22499559962325255,-0.10114691221770505,-0.03312123342077668,0.1452710030079996,-0.20877901085227965,0.20877901085227965,-0.14527100300799908,0.03312123342076048,0.1011469122175888,-0.22499559962349147,0.30648445130558705,-0.31960651385787014,0.24296094087234194,
                            0.1592678496914632,-0.2270309411372863,0.25664726816558986,-0.2561228730540197,0.2294466079023453,-0.1808325676183818,0.11550693724726287,-0.03969438377842774,-0.03969438377842774,0.11550693724726346,-0.18083256761839137,0.229446607902317,-0.2561228730539323,0.2566472681662554,-0.22703094113514646,0.15926784969576815,
                            -0.07789065669406012,0.11629958289316977,-0.14346081988963955,0.16405837863026962,-0.1796841248013294,0.1910402276970876,-0.1984721192331684,0.20215190531867153,-0.20215190531867153,0.19847211923316818,-0.19104022769707935,0.17968412480141777,-0.16405837862987688,0.14346081989075976,-0.11629958289061507,0.07789065669850084
                    }),
                    new double[16]{-0.9894009349916443,-0.9445750230731486,-0.865631202388051,-0.755404408354798,-0.617876244402738,-0.45801677765720394,-0.28160355077926175,-0.0950125098376374,0.09501250983763737,0.28160355077926164,0.4580167776571992,0.6178762444027687,0.7554044083547129,0.8656312023881924,0.944575023072992,0.9894009349917203},
                    new double[16]{0.02715245941180501,0.06225352393872681,0.09515851168238959,0.12462897125559844,0.14959598881654815,0.16915651939500606,0.18260341504492264,0.18945061045506853,0.1894506104550685,0.1826034150449232,0.1691565193950095,0.14959598881652297,0.12462897125576994,0.09515851168216734,0.06225352393886759,0.027152459411602402}
            },
            {
                    18,
                    array2d<double>(18, 18, new double[324]{
                            1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                            -0.9915651684212035,-0.9558239495710681,-0.8926024664973684,-0.8037049589730074,-0.6916870430600279,-0.5597708310740721,-0.4117511614628011,-0.2518862256915131,-0.08477501304173529,0.08477501304173529,0.25188622569151287,0.41175116146279866,0.5597708310741016,0.6916870430600094,0.8037049589728138,0.8926024664979532,0.9558239495703588,0.9915651684215105,
                            0.9748022248392544,0.8703991338604535,0.6951087447957786,0.46891249161670534,0.21764644830568738,-0.029984925017963945,-0.24569147155105145,-0.4048299939603262,-0.48921979574566044,-0.48921979574566044,-0.4048299939603264,-0.24569147155105445,-0.029984925017914353,0.21764644830564892,0.4689124916162385,0.6951087447973444,0.8703991338584196,0.9748022248401679,
                            -0.9499231081360295,-0.7493645966688471,-0.43902465581609057,-0.09230885207921977,0.2102193149382217,0.40115503137771374,0.44310702226295107,0.3378759825033164,0.12563936630217926,-0.12563936630217926,-0.3378759825033162,-0.4431070222629505,-0.40115503137768865,-0.21021931493826043,0.09230885207857219,0.4390246558187077,0.7493645966650508,0.9499231081378333,
                            0.9172419981062258,0.6006567494039686,0.1644488000137173,-0.2218534749094732,-0.4176952948307814,-0.3704823555182616,-0.1350186007076248,0.1546864599972435,0.34827545870480575,0.34827545870480575,0.15468645999724387,-0.13501860070762026,-0.3704823555182949,-0.4176952948307926,-0.22185347491006516,0.16444880001708023,0.6006567493982139,0.9172419981091813,
                            -0.8771709028945381,-0.43392811445819945,0.08700239654427468,0.39479560997362295,0.3518705101362842,0.052369363781955815,-0.2544162996215568,-0.3404348854363695,-0.15365659483861266,0.15365659483861266,0.34043488543636946,0.2544162996215603,-0.05236936378202918,-0.3518705101362532,-0.3947956099738839,-0.08700239654079213,0.4339281144505686,0.877170902898877,
                            0.8302138773599484,0.25984399651419104,-0.2794146818802009,-0.39683728502249144,-0.09812508758403048,0.2549914187389673,0.304568546507468,0.028304523710495,-0.2663481092297861,-0.2663481092297861,0.028304523710494526,0.30456854650746573,0.25499141873891695,-0.0981250875839699,-0.3968372850222425,-0.2794146818773976,0.2598439965050505,0.8302138773658669,
                            -0.7769599576176086,-0.08931111550251818,0.38860952365976215,0.25392108008167963,-0.17555556985209272,-0.3099705773838893,-0.014826584039141683,0.278560793885942,0.1736393152383475,-0.1736393152383475,-0.2785607938859422,0.014826584039135537,0.3099705773839138,0.17555556985214738,-0.2539210800809417,-0.3886095236584037,0.0893111154924913,0.7769599576252618,
                            0.7180749158693611,-0.06730280352451444,-0.40590056458558266,-0.035412934197997556,0.31353978851798703,0.10221842306094674,-0.25505085969658,-0.15632700886711084,0.20545407955002115,0.20545407955002115,-0.1563270088671104,-0.2550508596965828,0.1022184230610337,0.3135397885179988,-0.035412934197011185,-0.40590056458618795,-0.06730280353460519,0.7180749158788586,
                            -0.6542919569156543,0.2008991843571405,0.33892857527243503,-0.1719469196221416,-0.25359659973463317,0.1674492734885212,0.21154555149108112,-0.17323153411727465,-0.1872455389503781,0.1872455389503781,0.1732315341172751,-0.21154555149107668,-0.16744927348844535,0.253596599734589,0.17194691962299608,-0.3389285752751115,-0.20089918436635582,0.6542919569270562,
                            0.5864014931834232,-0.3042735553578485,-0.20949360815674317,0.2944413655410496,0.05109220651481082,-0.2700896968233484,0.06404793332231325,0.22360011884954648,-0.15474855988262143,-0.15474855988262143,0.2236001188495462,0.06404793332232024,-0.2700896968233554,0.051092206514733236,0.29444136554140343,-0.2094936081611141,-0.3042735553652317,0.5864014931967382,
                            -0.5152401484737861,0.3725890124094692,0.048872635068603204,-0.2954595001484746,0.163075439558714,0.13640566184337596,-0.24266032223096873,0.049959795558894585,0.19526822039370906,-0.19526822039370906,-0.0499597955588952,0.2426603222309699,-0.1364056618434676,-0.1630754395587781,0.2954595001481318,-0.04887263507385224,-0.3725890124141521,0.5152401484889696,
                            0.4416791517295152,-0.40366411862064366,0.10842345781463097,0.18523142369031673,-0.26302909575806555,0.10123339325929147,0.13279442770032918,-0.22908646226307927,0.11012460352558615,0.11012460352558615,-0.2290864622630793,0.13279442770032254,0.10123339325919183,-0.2630290957580736,0.18523142368935475,0.10842345780960273,-0.4036641186219483,0.4416791517464694,
                            -0.3666122907796033,0.39805597346838956,-0.23122675135469717,-0.013559334050519138,0.19934155095282033,-0.23488888146542772,0.11884364390282016,0.06485196628499786,-0.19820107786361124,0.19820107786361124,-0.0648519662849972,-0.1188436439028271,0.2348888814654108,-0.1993415509527625,0.013559334049279753,0.23122675135103316,-0.3980559734659143,0.3666122907981777,
                            0.2909436021089027,-0.35893536722330477,0.29736581428949715,-0.15098360710865893,-0.021673920675654632,0.15957374186458403,-0.21768184193997417,0.18121931785589956,-0.06985374099226656,-0.06985374099226656,0.18121931785589993,-0.21768184193997298,0.1595737418646717,-0.021673920675562886,-0.1509836071096918,0.2973658142881196,-0.358935367216986,0.2909436021288966,
                            -0.21557497613632925,0.2917671974204034,-0.2973519866901029,0.24725817437708325,-0.15706841202222074,0.0465351522388664,0.06236538477045061,-0.14877869186041678,0.1964365861597284,-0.1964365861597284,0.14877869186041626,-0.062365384770442144,-0.04653515223874665,0.15706841202229022,-0.24725817437747485,0.29735198669148155,-0.29176719741054474,0.21557497615749469,
                            0.1413938582179352,-0.20382436352835165,0.235465212782461,-0.2434779463167389,0.23081353498113272,-0.20007036088118316,0.15432362631600496,-0.09728471062507635,0.03322286100771784,0.03322286100771784,-0.09728471062507701,0.15432362631601088,-0.2000703608811381,0.2308135349811342,-0.24347794631628764,0.23546521278647364,-0.2038243635156172,0.1413938582399813,
                            -0.0692612236152485,0.10357598298631338,-0.12812962312875187,0.14714444083485004,-0.16208114985014263,0.17360145801422755,-0.18204487792708818,0.18759490957661404,-0.19034875519045202,0.19034875519045202,-0.18759490957661384,0.1820448779270842,-0.1736014580143028,0.16208114985007088,-0.14714444083368594,0.1281296231346743,-0.10357598297168366,0.06926122363784679
                    }),
                    new double[18]{-0.9915651684212035,-0.955823949571068,-0.8926024664973685,-0.8037049589730074,-0.6916870430600278,-0.5597708310740722,-0.411751161462801,-0.2518862256915131,-0.08477501304173532,0.08477501304173528,0.25188622569151287,0.41175116146279866,0.5597708310741016,0.6916870430600094,0.8037049589728138,0.8926024664979532,0.9558239495703588,0.9915651684215105},
                    new double[18]{0.021616013525702228,0.049714548896052264,0.07642573025453561,0.10094204410574469,0.12255520671150899,0.14064291467060788,0.15468467512626724,0.16427648374583237,0.1691423829631436,0.1691423829631436,0.16427648374583229,0.15468467512626863,0.1406429146706446,0.12255520671165103,0.10094204410661786,0.07642573025527827,0.04971454889663268,0.021616013524918303}
            }
    };

    template<class D>
    D *getCoefficients(int n, std::function<D(double)> V, double a, double b) {
        LegendreData &data = legendreData[(n + 1) / 2];

        double m = (a + b) / 2;
        double h = (b - a) / 2;
        D *aV = new D[data.n];
        for (int j = 0; j < data.n; ++j)
            aV[j] = V(m + data.nodes[j] * h);

        D *coeffs = new D[n];
        double H = 1;
        for (int i = 0; i < n; ++i) {
            coeffs[i] = aV[0] * data.weights[0] * data.polynomial[i][0] * (i + .5) / H;
            for (int j = 1; j < data.n; ++j)
                coeffs[i] += aV[j] * data.weights[j] * data.polynomial[i][j] * (i + .5) / H;
            H *= h * 2;
        }
		delete[] aV;

        return coeffs;
    }
}


#endif //SCHRODINGER_LEGENDRE_H
