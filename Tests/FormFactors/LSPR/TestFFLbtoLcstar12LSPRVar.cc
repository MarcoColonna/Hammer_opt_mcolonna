///
/// @file  TestFFLbtoLcstar12LSPRVar.cc
/// @brief Tests for FFLbtoLcstar12LSPRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                    \
    FRIEND_TEST(FFLbtoLcstar12LSPRVarTest, eval); \
    FRIEND_TEST(FFLbtoLcstar12LSPRVarTest, clone)

#include "Hammer/FormFactors/LSPR/FFLbtoLcstar12LSPRVar.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFLbtoLcstar12LSPRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3053.24065836299, 0, 0, -1608.80499684731};
        FourMomentum kNuTau{103.003318909325, -16.4797557212763, -28.5437742055744, -97.5876749994806};
        FourMomentum pTau{2463.75602272769, 16.4797557212763, 28.5437742055744, 1706.39267184679};

        // Evaluate the FF class
        FFLbtoLcstar12LSPRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -14122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({12, 7}, {FF_LBLCSTAR12, FF_LBLCSTAR12_VAR},
                                               {0.800689720050897,
                                                0.821582790732750,
                                                0.205904553691841,
                                                -0.138295143180212,
                                                -0.0244209855334506,
                                                -0.607291751416109,
                                                0.175868559662991,
                                                0.151804731061730,
                                                0.155765899659862,
                                                0.0390379501757416,
                                                -0.0761832938042431,
                                                -0.0134529027780799,
                                                -0.0492694608497764,
                                                -0.0142681818002912,
                                                0.164627941225841,
                                                0.168923716638095,
                                                0.0423355538536932,
                                                -0.107239218492228,
                                                -0.0189369441557653,
                                                -0.0492694608497764,
                                                0,
                                                -0.699710351366009,
                                                -0.717968482402067,
                                                -0.179936801989191,
                                                0.175868559662991,
                                                0.0310559246879845,
                                                0.558022290566333,
                                                -0.0808001889313500,
                                                0.00860587827688889,
                                                0.00883043870100303,
                                                0.00221307890104606,
                                                -0.175868559662991,
                                                -0.0310559246879845,
                                                0,
                                                -0.0808001889313500,
                                                0.811514188436004,
                                                0.832689711080709,
                                                0.208688162967660,
                                                -0.107239218492228,
                                                -0.0189369441557653,
                                                -0.607291751416109,
                                                0,
                                                -0.636956213012620,
                                                -0.653576847506181,
                                                -0.163799011623718,
                                                -0.0142681818002912,
                                                -0.00251956108740203,
                                                0.558022290566333,
                                                0.0808001889313500,
                                                0.00834708275883007,
                                                0.00856489021370339,
                                                0.00214652730895123,
                                                -0.0142681818002912,
                                                -0.00251956108740203,
                                                0,
                                                -0.0808001889313500,
                                                0.854156249543776,
                                                0.876444467374187,
                                                0.219653952013057,
                                                -0.0761832938042431,
                                                -0.0134529027780799,
                                                -0.607291751416109,
                                                -0.175868559662991,
                                                -0.687821574518799,
                                                -0.705769481695641,
                                                -0.176879496232211,
                                                -0.0142681818002912,
                                                -0.00251956108740203,
                                                0.558022290566333,
                                                0.0808001889313500,
                                                0.00839946552333471,
                                                0.00861863984576471,
                                                0.00215999800737082,
                                                -0.0142681818002912,
                                                -0.00251956108740203,
                                                0,
                                                0.0808001889313500,
                                                -0.0300044513272222,
                                                -0.0307873827258044,
                                                -0.00771591417322995,
                                                0.190136741463282,
                                                0.0335754857753865,
                                                0,
                                                -0.161600377862700})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 7; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFLbtoLcstar12LSPRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcstar12LSPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3053.24065836299, 0, 0, -1608.80499684731};
        FourMomentum kNuTau{103.003318909325, -16.4797557212763, -28.5437742055744, -97.5876749994806};
        FourMomentum pTau{2463.75602272769, 16.4797557212763, 28.5437742055744, 1706.39267184679};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pLbmes, -5122), {Particle(pLcmes, -14122), Particle(kNuTau, 16), Particle(pTau, -15)},
                     {});

        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -14122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 7; ++idx2) {
                EXPECT_NEAR(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real(), 1e-12);
            }
        }
    }


} // namespace Hammer
