///
/// @file  TestFFBtoDBCL.cc
/// @brief Tests for FFBtoDBCL
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                 \
    FRIEND_TEST(FFBtoDBCLTest, eval);          \
    FRIEND_TEST(FFBtoDBCLTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDBCLTest, uninitialized); \
    FRIEND_TEST(FFBtoDBCLTest, clone);         \
    FRIEND_TEST(FFBtoDBCLTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BCL/FFBtoDBCL.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoDBCLTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 =  0.8958 GeV^2)
        FourMomentum pBmes{8152.95906902087, 0, 0, 6212.27346316384};
        FourMomentum pDmes{1869., 0, 0, 0};
        FourMomentum kNuMu{102.442823922292, 30.030225436543, 92.423530443711, 32.4130954048573};
        FourMomentum pMu{6181.51624509857, -30.030225436543, -92.423530443711, 6179.86036775898};


        // Evaluate the FF class
        FFBtoDBCL ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        // ff.setSignatureIndex(1); //Adjust to B to D- signature
        ff.initSettings();
        set.changeSetting<bool>("BtoDBCL", "WithFpPole", true);
        set.changeSetting<bool>("BtoDBCL", "WithF0Pole", true);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 521), {Particle(pDmes, -421), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BD}, {0., 0.7140215390817, 0.7274761053638, 0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 3; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoDBCLTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBCL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BPLUS);
        const double Mc = pdg.getMass(PID::D0);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoDBCLTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBCL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BPLUS);
        const double Mc = pdg.getMass(PID::D0);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDBCLTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBCL ff;
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
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFBtoDBCLTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBCL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Belle-II:2025rna"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Belle-II:2025rna"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Belle-II:2025rna"));

        const string tmpfile = testing::TempDir() + "hammer_bcl_btod_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Belle-II:2025rna"), string::npos);
    }

} // namespace Hammer
