///
/// @file  TestFFBtoD0starBLR.cc
/// @brief Tests for FFBtoD0starBLR
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                      \
    FRIEND_TEST(FFBtoD0starBLRTest, eval);          \
    FRIEND_TEST(FFBtoD0starBLRTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoD0starBLRTest, uninitialized); \
    FRIEND_TEST(FFBtoD0starBLRTest, clone);         \
    FRIEND_TEST(FFBtoD0starBLRTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLR/FFBtoD0starBLR.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/SettingsConsumer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoD0starBLRTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD0starmes{3117.27272727273, 0, 0, -2104.13622567512};
        FourMomentum kNuTau{49.97993502618, -24.8428481565870, -43.029075211928, -5.41530120089};
        FourMomentum pTau{2112.74733770109, 24.8428481565870, 43.029075211928, 2109.55152687601};

        // Evaluate the FF class
        FFBtoD0starBLR ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pD0starmes, -10411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval",
                      MD::makeVector({4}, {FF_BDSSD0STAR},
                                     {0.495949833064184, -0.131663927760578, 0.720459193477596, 0.877081840967645})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoD0starBLRTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD0starBLR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD0STARPLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoD0starBLRTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD0starBLR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD0STARPLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD0starBLRTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starBLR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD0starmes{3117.27272727273, 0, 0, -2104.13622567512};
        FourMomentum kNuTau{49.97993502618, -24.8428481565870, -43.029075211928, -5.41530120089};
        FourMomentum pTau{2112.74733770109, 24.8428481565870, 43.029075211928, 2109.55152687601};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pD0starmes, -10411), Particle(kNuTau, 14), Particle(pTau, -13)},
                     {});

        ff.eval(Particle(pBmes, 511), {Particle(pD0starmes, -10411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFBtoD0starBLRTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starBLR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Bernlochner:2016bci"));
        EXPECT_FALSE(set.checkReference("Bernlochner:2017jxt"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2016bci"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2017jxt"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2016bci"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2017jxt"));

        const string tmpfile = testing::TempDir() + "hammer_blr_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Bernlochner:2016bci"), string::npos);
        EXPECT_NE(content.find("Bernlochner:2017jxt"), string::npos);
    }

} // namespace Hammer
