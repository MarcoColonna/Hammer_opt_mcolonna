///
/// @file  TestFFBtoDstarBGLVar.cc
/// @brief Tests for FFBtoDstarBGLVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS               \
    FRIEND_TEST(FFBtoDstarBGLVarTest, eval); \
    FRIEND_TEST(FFBtoDstarBGLVarTest, clone)

#include "Hammer/FormFactors/BGL/FFBtoDstarBGLVar.hh"
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

    TEST(FFBtoDstarBGLVarTest, eval) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBGLVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({8, 11}, {FF_BDSTAR, FF_BDSTAR_VAR},
                                               {0.,
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
                                                5057.250417739538,
                                                0.,
                                                0.,
                                                0.,
                                                9.74419145803456e6,
                                                148284.18467327964,
                                                2256.5442724433515,
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
                                                0.00022750456342981215,
                                                0.,
                                                0.,
                                                0.,
                                                -0.6155933683815917,
                                                0.03948027970144063,
                                                0.000600799061822005,
                                                -0.2912197931470863,
                                                -0.004431695515581247,
                                                0.08433248804219676,
                                                0.0012833465233788206,
                                                -0.00010514774368697381,
                                                0.,
                                                0.,
                                                0.,
                                                -0.2007283073085566,
                                                -0.01954599965590657,
                                                -0.0002974451636221157,
                                                0.09831714830739159,
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
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 5; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBtoDstarBGLVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarBGLVar ff;
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
