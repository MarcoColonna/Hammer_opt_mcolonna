///
/// @file  TestFFBstoKBCLVar.cc
/// @brief Tests for FFBstoKBCLVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS              \
    FRIEND_TEST(FFBstoKBCLVarTest, eval);   \
    FRIEND_TEST(FFBstoKBCLVarTest, clone);  \
    FRIEND_TEST(FFBstoKBCLVarTest, addRefs)

#include "Hammer/FormFactors/BCL/FFBstoKBCLVar.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBstoKBCLVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 = 19.8958 GeV^2)
        FourMomentum pBsmes{9269.42342515698, 0, 0, 7557.61348805621};
        FourMomentum pKmes{493.7, 0, 0, 0};
        FourMomentum kNuMu{691.247703926747, 143.20931144041, 440.75294019035, -512.885296157276};
        FourMomentum pMu{8084.47572123023, -143.20931144041, -440.75294019035, 8070.49878421349};

        // Evaluate the FF class
        FFBstoKBCLVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        double sqtp = 5.280 + 0.139; // GeV t+ = (mB + mPi)^2 branch point
        set.changeSetting<double>("BstoKBCLVar", "tp", sqtp * sqtp);
        ff.calcUnits();
        ff.eval(Particle(pBsmes, 531), {Particle(pKmes, -321), Particle(kNuMu, 14), Particle(pMu, -13)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({4, 8}, {FF_BSK, FF_BSK_VAR},
                                               {0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.5601214962546,
                                                -0.0010217667629501596,
                                                0.0006708994336375612,
                                                0.001042098943669952,
                                                0.0018240266408726116,
                                                -0.007149561500480032,
                                                0.01096659699170062,
                                                0.00083984895240395,
                                                1.4242702517186,
                                                0.007733427260543903,
                                                0.006438654669961988,
                                                0.011127399375760485,
                                                -0.016500091853727375,
                                                -0.028524289922883937,
                                                0.011418470166546246,
                                                -0.0009615519864016019,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBstoKBCLVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBstoKBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBsmes{9269.42342515698, 0, 0, 7557.61348805621};
        FourMomentum pKmes{493.7, 0, 0, 0};
        FourMomentum kNuMu{691.247703926747, 143.20931144041, 440.75294019035, -512.885296157276};
        FourMomentum pMu{8084.47572123023, -143.20931144041, -440.75294019035, 8070.49878421349};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBsmes, 531), {Particle(pKmes, -321), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBsmes, 531), {Particle(pKmes, -321), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBstoKBCLVarTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBstoKBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("FlavourLatticeAveragingGroupFLAG:2021npn"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("FlavourLatticeAveragingGroupFLAG:2021npn"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("FlavourLatticeAveragingGroupFLAG:2021npn"));

        const string tmpfile = testing::TempDir() + "hammer_bcl_bstokvar_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("FlavourLatticeAveragingGroupFLAG:2021npn"), string::npos);
    }

} // namespace Hammer
