///
/// @file  TestFFLbtoLcBLRS.cc
/// @brief Tests for FFLbtoLcBLRS
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                    \
    FRIEND_TEST(FFLbtoLcBLRSTest, eval);          \
    FRIEND_TEST(FFLbtoLcBLRSTest, zeroRecoil);    \
    FRIEND_TEST(FFLbtoLcBLRSTest, uninitialized); \
    FRIEND_TEST(FFLbtoLcBLRSTest, clone);         \
    FRIEND_TEST(FFLbtoLcBLRSTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLRS/FFLbtoLcBLRS.hh"
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

    TEST(FFLbtoLcBLRSTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{2514.6000000000004, 380.2511160960518, -878.708483534587, -425.0854616098454};
        FourMomentum pTau{2456.13415321897, -581.6942711749232, 1418.910848905593, 723.6624190525239};
        FourMomentum kNuTau{649.2658467810297, 201.44315507887137, -540.2023653710064, -298.5769574426786};

        // Evaluate the FF class
        FFLbtoLcBLRS ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({12}, {FF_LBLC},
                                               {0.7312555660768754, 1.1660722177550018, 1.195461806802338,
                                                -0.31472386413207737, -0.0851912009508891, 0.7527695959975337,
                                                -0.3856860143396382, 0.11604815859155679, 0.805646707782559,
                                                -0.37394516263654165, 0.11780903740914717, 0.013524178395048072})};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFLbtoLcBLRSTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcBLRS ff;
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

    TEST(FFLbtoLcBLRSTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcBLRS ff;
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

    TEST(FFLbtoLcBLRSTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcBLRS ff;
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

    TEST(FFLbtoLcBLRSTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcBLRS ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Bernlochner:2018kxh"));
        EXPECT_FALSE(set.checkReference("Bernlochner:2018bfn"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2018kxh"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2018bfn"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2018kxh"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2018bfn"));

        const string tmpfile = testing::TempDir() + "hammer_blrs_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Bernlochner:2018kxh"), string::npos);
        EXPECT_NE(content.find("Bernlochner:2018bfn"), string::npos);
    }


} // namespace Hammer
