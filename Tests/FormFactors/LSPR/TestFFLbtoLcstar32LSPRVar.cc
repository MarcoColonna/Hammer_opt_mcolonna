///
/// @file  TestFFLbtoLcstar32LSPRVar.cc
/// @brief Tests for FFLbtoLcstar32LSPRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                    \
    FRIEND_TEST(FFLbtoLcstar32LSPRVarTest, eval); \
    FRIEND_TEST(FFLbtoLcstar32LSPRVarTest, clone)

#include "Hammer/FormFactors/LSPR/FFLbtoLcstar32LSPRVar.hh"
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

    TEST(FFLbtoLcstar32LSPRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3067.17304270463, 0, 0, -1586.48210638947};
        FourMomentum kNuTau{103.858273281499, -16.4797557212763, -28.5437742055744, -98.4896519157352};
        FourMomentum pTau{2448.96868401388, 16.4797557212763, 28.5437742055744, 1684.97175830520};

        // Evaluate the FF class
        FFLbtoLcstar32LSPRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4124), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({16, 7}, {FF_LBLCSTAR32, FF_LBLCSTAR32_VAR},
                                               {0.720401386256984,
                                                0.739199425878894,
                                                0.173078896728211,
                                                0.0535928980610493,
                                                0.00902755611543353,
                                                0.246712412519125,
                                                -0.0714466424738093,
                                                1.24417777068813,
                                                1.27664314829049,
                                                0.298917964335578,
                                                0.380054511070786,
                                                0.0640189941157290,
                                                0.246712412519125,
                                                0.0714466424738093,
                                                1.26182697446594,
                                                1.29475288759507,
                                                0.303158245901200,
                                                0.380054511070786,
                                                0.0640189941157290,
                                                0.246712412519125,
                                                0.0714466424738093,
                                                -0.470732504336548,
                                                -0.483015723714866,
                                                -0.113095094011401,
                                                -0.493424825038250,
                                                -0.0831158690412042,
                                                0,
                                                0,
                                                -0.0170396347149046,
                                                -0.0174842642431428,
                                                -0.00409383051361224,
                                                0.166963212028513,
                                                0.0281244310409087,
                                                0,
                                                -0.142893284947619,
                                                0.0207258015212351,
                                                0.0212666173020268,
                                                0.00497944469504894,
                                                -0.0521943581218030,
                                                -0.00879197643532670,
                                                0,
                                                0.0240699270808944,
                                                0.738050590034792,
                                                0.757309165183478,
                                                0.177319178293833,
                                                0.0535928980610493,
                                                0.00902755611543353,
                                                0.246712412519125,
                                                -0.0714466424738093,
                                                -0.511494766669524,
                                                -0.524841630062182,
                                                -0.122888366938587,
                                                -0.493424825038250,
                                                -0.0831158690412042,
                                                0,
                                                0,
                                                0.0281547899021519,
                                                0.0288894565285962,
                                                0.00676428455974778,
                                                -0.166963212028513,
                                                -0.0281244310409087,
                                                0,
                                                0.142893284947619,
                                                -0.0207258015212351,
                                                -0.0212666173020268,
                                                -0.00497944469504894,
                                                0.0521943581218030,
                                                0.00879197643532670,
                                                0,
                                                -0.309856496976132,
                                                0.781633411513197,
                                                0.802029229899632,
                                                0.187790100201657,
                                                0.0535928980610493,
                                                0.00902755611543353,
                                                0.246712412519125,
                                                -0.0714466424738093,
                                                -0.496146403812827,
                                                -0.509092769457124,
                                                -0.119200870272841,
                                                -0.493424825038250,
                                                -0.0831158690412042,
                                                0,
                                                0,
                                                0.0276299806183191,
                                                0.0283509529544696,
                                                0.00663819733452677,
                                                -0.166963212028513,
                                                -0.0281244310409087,
                                                0,
                                                0.142893284947619,
                                                -0.0207258015212351,
                                                -0.0212666173020268,
                                                -0.00497944469504894,
                                                0.0521943581218030,
                                                0.00879197643532670,
                                                0,
                                                -0.309856496976132,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0.285786569895237,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 16; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 7; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFLbtoLcstar32LSPRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcstar32LSPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3067.17304270463, 0, 0, -1586.48210638947};
        FourMomentum kNuTau{103.858273281499, -16.4797557212763, -28.5437742055744, -98.4896519157352};
        FourMomentum pTau{2448.96868401388, 16.4797557212763, 28.5437742055744, 1684.97175830520};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4124), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4124), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 16; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 7; ++idx2) {
                EXPECT_NEAR(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real(), 1e-12);
            }
        }
    }


} // namespace Hammer
