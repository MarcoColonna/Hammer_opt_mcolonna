///
/// @file  TestFFLbtoLcBLRSXP.cc
/// @brief Tests for FFLbtoLcBLRSXP
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                      \
    FRIEND_TEST(FFLbtoLcBLRSXPTest, eval);          \
    FRIEND_TEST(FFLbtoLcBLRSXPTest, zeroRecoil);    \
    FRIEND_TEST(FFLbtoLcBLRSXPTest, uninitialized); \
    FRIEND_TEST(FFLbtoLcBLRSXPTest, clone);         \
    FRIEND_TEST(FFLbtoLcBLRSXPTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLRSXP/FFLbtoLcBLRSXP.hh"
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

    TEST(FFLbtoLcBLRSXPTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{2514.6000000000004, 380.2511160960518, -878.708483534587, -425.0854616098454};
        FourMomentum pTau{2456.13415321897, -581.6942711749232, 1418.910848905593, 723.6624190525239};
        FourMomentum kNuTau{649.2658467810297, 201.44315507887137, -540.2023653710064, -298.5769574426786};

        // Evaluate the FF class
        FFLbtoLcBLRSXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({12}, {FF_LBLC},
                                               {0.744050262639, 1.18511593540, 1.21407897674, -0.301188694543,
                                                -0.0826568041477, 0.763812946280, -0.364472049319, 0.104855753419,
                                                0.830515836674, -0.352357783478, 0.1066387974459, 0.00453549559603})};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            double comRe = compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 2e-4);
        }
    }

    TEST(FFLbtoLcBLRSXPTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcBLRSXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACPLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 12; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFLbtoLcBLRSXPTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcBLRSXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACPLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFLbtoLcBLRSXPTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcBLRSXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{2514.6000000000004, 380.2511160960518, -878.708483534587, -425.0854616098454};
        FourMomentum kNuTau{649.2658467810297, 201.44315507887137, -540.2023653710064, -298.5769574426786};
        FourMomentum pTau{2456.13415321897, -581.6942711749232, 1418.910848905593, 723.6624190525239};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 12; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFLbtoLcBLRSXPTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcBLRSXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Bernlochner:2022ywh"));
        EXPECT_FALSE(set.checkReference("Bernlochner:2023jkp"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2022ywh"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2023jkp"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2022ywh"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2023jkp"));

        const string tmpfile = testing::TempDir() + "hammer_blrsxp_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Bernlochner:2022ywh"), string::npos);
        EXPECT_NE(content.find("Bernlochner:2023jkp"), string::npos);
    }


} // namespace Hammer
