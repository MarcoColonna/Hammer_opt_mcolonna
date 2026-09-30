///
/// @file  TestFFBtoOmegaBSZVar.cc
/// @brief Tests for FFBtoOmegaBSZVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS               \
    FRIEND_TEST(FFBtoOmegaBSZVarTest, eval); \
    FRIEND_TEST(FFBtoOmegaBSZVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BSZ/FFBtoOmegaBSZVar.hh"
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

    TEST(FFBtoOmegaBSZVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{17576.6777493606, 0, 0, 16764.8799788388};
        FourMomentum pOmegames{782., 0, 0, 0};
        FourMomentum kNuMu{195.8915304338, 31.769954441031, 97.777865759968, 166.7442672416};
        FourMomentum pMu{16598.7862189268, -31.769954441031, -97.777865759968, 16598.1357115972};

        // Evaluate the FF class
        FFBtoOmegaBSZVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 521), {Particle(pOmegames, 223), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Ap,V,A0,A1,A12,T1,T2,T23
        Tensor ffEval{
            "ffEval",
            MD::makeVector(
                {8, 9}, {FF_BOMEGA, FF_BOMEGA_VAR},
                {-0.11262003760786,      -0.005834456604024177,  0.002679301333039777,   0.0035775573146732144,
                 0.0006623063208313079,  0.002468493753665233,   0.007099732121495844,   -0.0005302335677653433,
                 -0.00260627977833401,   0.32133945171369,       -0.002320560975399416,  -0.016392110913949716,
                 0.018724801272900236,   -0.002715443282962616,  0.006821575545656255,   0.002243575749707538,
                 0.0049602298301523445,  0.0017069161345335255,  0.34599698254843,       0.01792491303218425,
                 -0.008231485233539724,  -0.01099115274744552,   -0.0020347710176435,    -0.0075838315130188245,
                 -0.021812156549736263,  0.001629010417591096,   0.008007144715405957,   0.24779174538358,
                 0.00019984469742809213, -0.013079576428950978,  0.014374390618791356,   -0.004043557736467908,
                 0.005206797741210745,   0.001999387707991577,   0.00021709066477626582, 0.0023382944777674162,
                 0.2743143152703,        0.014790505976363768,   -0.006707607941744271,  -0.008926185795994442,
                 -0.0017217723460563406, -0.006039541820020712,  -0.017759532095722584,  0.001299711400072675,
                 0.006359452109469523,   0.26569423684374,       -0.0024813241502187255, -0.015563714842465144,
                 0.01284925159123104,    -0.0024923472637945034, 0.004925454598639538,   0.0025279624977015804,
                 -0.003854287101873415,  0.0004114515455979505,  0.25639515670707,       -0.0023520986286957495,
                 -0.015564418958370883,  0.01275089009999633,    -0.0025149500656942714, 0.00482992872793258,
                 0.002576304930965115,   -0.00385806924580965,   0.0003965513190370141,  0.69258666245592,
                 0.032195179251141216,   -0.015142214228246714,  -0.020105302410552143,  -0.010987695517747163,
                 0.00036444798928915125, -0.04639689662873461,   -0.0018770420318346597, 0.0110055529425272})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; idx1++) {
            for (IndexType idx2 = 0; idx2 < 9; idx2++) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-4);
            }
        }
    }

    TEST(FFBtoOmegaBSZVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoOmegaBSZVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{17576.6777493606, 0, 0, 16764.8799788388};
        FourMomentum pOmegames{782., 0, 0, 0};
        FourMomentum kNuMu{195.8915304338, 31.769954441031, 97.777865759968, 166.7442672416};
        FourMomentum pMu{16598.7862189268, -31.769954441031, -97.777865759968, 16598.1357115972};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 521), {Particle(pOmegames, 223), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBmes, 521), {Particle(pOmegames, 223), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

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
