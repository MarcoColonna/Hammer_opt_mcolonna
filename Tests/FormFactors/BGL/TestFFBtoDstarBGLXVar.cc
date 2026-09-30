///
/// @file  TestFFBtoDstarBGLXVar.cc
/// @brief Tests for FFBtoDstarBGLXVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS               \
    FRIEND_TEST(FFBtoDstarBGLXVarTest, eval); \
    FRIEND_TEST(FFBtoDstarBGLXVarTest, evalOpt); \
    FRIEND_TEST(FFBtoDstarBGLXVarTest, clone)

#include "Hammer/FormFactors/BGL/FFBtoDstarBGLXVar.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoDstarBGLXVarTest, eval) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBGLXVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({8, 20}, {FF_BDSTAR, FF_BDSTAR_VAR},
                                               {-0.6845787512332282,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -106.10066319832539,
                                                -1.6146080876402802,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                5057.250417739539,
                                                0.,
                                                0.,
                                                0.,
                                                9.744191458034564e6,
                                                148284.18467327973,
                                                2256.544272443353,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.00011113797765262371,
                                                0.14078216488434997,
                                                0.002142380784113388,
                                                0.000032602108568999,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.0002275045634298128,
                                                0.,
                                                0.,
                                                0.,
                                                -0.6155933683815907,
                                                0.03948027970144065,
                                                0.0006007990618220053,
                                                -0.29121979314708624,
                                                -0.004431695515581247,
                                                0.08433248804219673,
                                                0.0012833465233788206,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -0.00010514774368697405,
                                                0.,
                                                0.,
                                                0.,
                                                -0.2007283073085571,
                                                -0.01954599965590658,
                                                -0.0002974451636221159,
                                                0.09831714830739158,
                                                0.0014961608912294603,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                7.948449742654126e-11,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -0.0002710222553379486,
                                                -4.124335439649424e-6,
                                                -6.276290040291185e-8,
                                                0.002419249579481028,
                                                0.00003681541490225166,
                                                5.60245948132266e-7,
                                                -0.0005000034982152994,
                                                -7.608903353957512e-6,
                                                -1.1579001038296013e-7,
                                                0.4019876405768768,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                6463.975918621697,
                                                98.3668478773781,
                                                1.4969172043874892,
                                                -18353.44611326193,
                                                -279.29724129198917,
                                                -4.250261695374419,
                                                0.,
                                                0.,
                                                0.,
                                                -0.9109755230439129,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -2182.2681493546575,
                                                -33.20910253653748,
                                                -0.5053661675850396,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 20; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBtoDstarBGLXVarTest, evalOpt) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBGLXVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<bool>("BtoD*BGLXVar", "OptZ", true);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({8, 20}, {FF_BDSTAR, FF_BDSTAR_VAR},
                                               {-0.7921776782898919,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -106.10405848178479,
                                                1.3735908032610662,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                5616.443068618494,
                                                0.,
                                                0.,
                                                0.,
                                                9.744503277877023e6,
                                                -126149.36955627412,
                                                1633.0912911255587,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                4.462384142511799e-6,
                                                0.14078666999620687,
                                                -0.0018225813215404876,
                                                0.000023594582311789653,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.00036546970467231854,
                                                0.,
                                                0.,
                                                0.,
                                                -0.6156130677304741,
                                                0.05680149051595911,
                                                -0.0021103374614617684,
                                                -0.2911229526643424,
                                                0.011966184945240366,
                                                0.0843351867312822,
                                                -0.001091777624182782,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -0.00014629821908523222,
                                                0.,
                                                0.,
                                                0.,
                                                -0.20073473073213544,
                                                -0.013887241857390337,
                                                0.0006439872601709775,
                                                0.09828445451277971,
                                                -0.004039839350276132,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -4.234181740041228e-8,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -0.0002710309282090973,
                                                3.5086837931956938e-6,
                                                -4.5422351028279814e-8,
                                                0.0024193269968866607,
                                                -0.0000313198699517722,
                                                4.054574908882719e-7,
                                                -0.0005000194986205346,
                                                6.473100035794435e-6,
                                                -8.379878022557007e-8,
                                                1.216137728315007,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                6464.1827696427745,
                                                -83.68334001572875,
                                                1.083339015269685,
                                                -18354.03343430344,
                                                237.60572299346467,
                                                -3.075971273634638,
                                                0.,
                                                0.,
                                                0.,
                                                -1.0887944086920733,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                -2182.337983215522,
                                                28.251882409066933,
                                                -0.3657402592057256,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 20; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }
    
    TEST(FFBtoDstarBGLXVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarBGLXVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }


} // namespace Hammer
