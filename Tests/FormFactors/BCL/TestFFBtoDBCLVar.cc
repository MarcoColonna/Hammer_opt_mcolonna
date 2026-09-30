///
/// @file  TestFFBtoDBCLVar.cc
/// @brief Tests for FFBtoDBCLVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS             \
    FRIEND_TEST(FFBtoDBCLVarTest, eval);   \
    FRIEND_TEST(FFBtoDBCLVarTest, clone);  \
    FRIEND_TEST(FFBtoDBCLVarTest, addRefs)

#include "Hammer/FormFactors/BCL/FFBtoDBCLVar.hh"
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

    TEST(FFBtoDBCLVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 =  0.8958 GeV^2)
        FourMomentum pBmes{8152.95906902087, 0, 0, 6212.27346316384};
        FourMomentum pDmes{1869., 0, 0, 0};
        FourMomentum kNuMu{102.442823922292, 30.030225436543, 92.423530443711, 32.4130954048573};
        FourMomentum pMu{6181.51624509857, -30.030225436543, -92.423530443711, 6179.86036775898};

        // Evaluate the FF class
        FFBtoDBCLVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        // ff.setSignatureIndex(1); //Adjust to B to D- signature
        ff.initSettings();
        set.changeSetting<bool>("BtoDBCLVar", "WithFpPole", true);
        set.changeSetting<bool>("BtoDBCLVar", "WithF0Pole", true);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 521), {Particle(pDmes, -421), Particle(kNuMu, 14), Particle(pMu, -13)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({4, 6}, {FF_BD, FF_BD_VAR},
                                               {0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.7140215390817,
                                                0.000013156563695557199,
                                                0.008253770267909338,
                                                -0.0029739090631140082,
                                                0.005917869980104652,
                                                0.0003614103429968246,
                                                0.7274761053638,
                                                0.00033592519216466346,
                                                0.008501208972565474,
                                                -0.0032958261978939933,
                                                0.00599301259201318,
                                                0.0007560357283222125,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBtoDBCLVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{8152.95906902087, 0, 0, 6212.27346316384};
        FourMomentum pDmes{1869., 0, 0, 0};
        FourMomentum kNuMu{102.442823922292, 30.030225436543, 92.423530443711, 32.4130954048573};
        FourMomentum pMu{6181.51624509857, -30.030225436543, -92.423530443711, 6179.86036775898};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 521), {Particle(pDmes, -421), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBmes, 521), {Particle(pDmes, -421), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBtoDBCLVarTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBCLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Belle-II:2025rna"));
        EXPECT_FALSE(set.checkReference("FlavourLatticeAveragingGroupFLAG:2024oxs"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Belle-II:2025rna"));
        EXPECT_TRUE(set.checkReference("FlavourLatticeAveragingGroupFLAG:2024oxs"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Belle-II:2025rna"));
        EXPECT_TRUE(set.checkReference("FlavourLatticeAveragingGroupFLAG:2024oxs"));

        const string tmpfile = testing::TempDir() + "hammer_bcl_btodvar_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Belle-II:2025rna"), string::npos);
        EXPECT_NE(content.find("FlavourLatticeAveragingGroupFLAG:2024oxs"), string::npos);
    }

} // namespace Hammer
