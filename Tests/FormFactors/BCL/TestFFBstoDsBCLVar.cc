///
/// @file  TestFFBstoDsBCLVar.cc
/// @brief Tests for FFBstoDsBCLVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS               \
    FRIEND_TEST(FFBstoDsBCLVarTest, eval);   \
    FRIEND_TEST(FFBstoDsBCLVarTest, clone);  \
    FRIEND_TEST(FFBstoDsBCLVarTest, addRefs)

#include "Hammer/FormFactors/BCL/FFBstoDsBCLVar.hh"
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

    TEST(FFBstoDsBCLVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (w =1.1)
        FourMomentum pBsmes{5903.59, 0, 0, 2459.42254972584};
        FourMomentum pDsmes{1968.3, 0, 0, 0};
        FourMomentum kNuMu{763.91243505469, 98.5729691755825, 303.376404442174, -694.124149851331};
        FourMomentum pMu{3171.37756494531, -98.5729691755825, -303.376404442174, 3153.54669957717};

        // Evaluate the FF class
        FFBstoDsBCLVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBsmes, 531), {Particle(pDsmes, -431), Particle(kNuMu, 14), Particle(pMu, -13)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({4, 7}, {FF_BSDS, FF_BSDS_VAR},
                                               {0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.8582085431933,
                                                -8.50533670891281e-6,
                                                0.00023142056890030314,
                                                0.00025571054862531817,
                                                -0.00385014787695205,
                                                0.01172226072152934,
                                                0.,
                                                1.0749552659518,
                                                -0.002620554774529168,
                                                -0.00634058349760404,
                                                0.026743190845257794,
                                                0.011006098918111096,
                                                0.012074460125540773,
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
            for (IndexType idx2 = 0; idx2 < 7; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBstoDsBCLVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBstoDsBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBsmes{5903.59, 0, 0, 2459.42254972584};
        FourMomentum pDsmes{1968.3, 0, 0, 0};
        FourMomentum kNuMu{763.91243505469, 98.5729691755825, 303.376404442174, -694.124149851331};
        FourMomentum pMu{3171.37756494531, -98.5729691755825, -303.376404442174, 3153.54669957717};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBsmes, 531), {Particle(pDsmes, -431), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBsmes, 531), {Particle(pDsmes, -431), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 7; ++idx2) {
                EXPECT_NEAR(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real(), 1e-12);
            }
        }
    }

    TEST(FFBstoDsBCLVarTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBstoDsBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("McLean:2019qcx"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("McLean:2019qcx"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("McLean:2019qcx"));

        const string tmpfile = testing::TempDir() + "hammer_bcl_bstodsvar_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("McLean:2019qcx"), string::npos);
    }

} // namespace Hammer
