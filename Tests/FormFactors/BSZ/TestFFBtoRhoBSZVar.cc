///
/// @file  TestFFBtoRhoBSZVar.cc
/// @brief Tests for FFBtoRhoBSZVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS             \
    FRIEND_TEST(FFBtoRhoBSZVarTest, eval); \
    FRIEND_TEST(FFBtoRhoBSZVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BSZ/FFBtoRhoBSZVar.hh"
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

    TEST(FFBtoRhoBSZVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{17838.5064935065, 0, 0, 17039.1875956242};
        FourMomentum pRhomes{770., 0, 0, 0};
        FourMomentum kNuMu{198.6188302152, 31.769954441031, 97.777865759968, 169.9399855193};
        FourMomentum pMu{16869.8876632913, -31.769954441031, -97.777865759968, 16869.2476101049};

        // Evaluate the FF class
        FFBtoRhoBSZVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 521), {Particle(pRhomes, 113), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Ap,V,A0,A1,A12,T1,T2,T23
        Tensor ffEval{
            "ffEval",
            MD::makeVector(
                {8, 9}, {FF_BRHO, FF_BRHO_VAR},
                {-0.12031819247473,       -0.0016459197011122423,  0.00350691292444776,    -0.002777250825231299,
                 -0.0018840892990992888,  0.0018400219786380721,   -0.007914206538683563,  0.0016136674842608801,
                 0.0008889990468018276,   0.34558307298139,        -0.003054217607980781,  -0.0024592748298704146,
                 -0.015566915368680252,   -0.005810348105715765,   0.010573973105749713,   0.0010642622808846593,
                 0.004371554555962617,    -0.00017575520409936104, 0.37540838626044,       0.005135483223275534,
                 -0.01094202376751395,    0.008665383256647008,    0.0058786033001117415,  -0.0057411075372441125,
                 0.024693352219723713,    -0.0050348521180996935,  -0.0027737924804433645, 0.26721909723634,
                 -0.00010560779639919283, -0.002260520064704657,   -0.013430259314281684,  -0.0055688028557477134,
                 0.0076553105395672245,   -0.0006851606813036665,  -0.0012832746789069372, -0.0032600845229213923,
                 0.30077407846314,        0.004380661930302681,    -0.008850412743044885,  0.007232294209802881,
                 0.0046298780757009295,   -0.00461152597664372,    0.020290209184163542,   -0.004160152113532799,
                 -0.002035927125815942,   0.2871289672528,         -0.0055917516873967826, -0.004870911704324534,
                 -0.01190622872494191,    -0.004549549853839193,   0.008158295919339146,   -0.0018246295639384887,
                 -0.004762790734380338,   -0.0007908649506400517,  0.2770965126059,        -0.0055252966866519825,
                 -0.004985180005887252,   -0.011909414058292198,   -0.004542599803578778,  0.00807957639219887,
                 -0.0018878263456153699,  -0.004780593035661231,   -0.0007688003768989505, 0.7574294922472,
                 0.008891758503701847,    -0.007588005190628789,   0.015215738577649424,   -0.003939921200847206,
                 -0.00044098229876847125, 0.04698887228154348,     -0.011617780861640833,  0.004260362184917721})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; idx1++) {
            for (IndexType idx2 = 0; idx2 < 9; idx2++) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-4);
            }
        }
    }

    TEST(FFBtoRhoBSZVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoRhoBSZVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{17838.5064935065, 0, 0, 17039.1875956242};
        FourMomentum pRhomes{770., 0, 0, 0};
        FourMomentum kNuMu{198.6188302152, 31.769954441031, 97.777865759968, 169.9399855193};
        FourMomentum pMu{16869.8876632913, -31.769954441031, -97.777865759968, 16869.2476101049};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 521), {Particle(pRhomes, 113), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBmes, 521), {Particle(pRhomes, 113), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 9; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

} // namespace Hammer
