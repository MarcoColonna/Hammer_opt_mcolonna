///
/// @file  TestFFBtoD0starLLSW.cc
/// @brief Unit tests for the FFBtoD0starLLSW form factor class (B -> D0*)
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                       \
    FRIEND_TEST(FFBtoD0starLLSWTest, eval);          \
    FRIEND_TEST(FFBtoD0starLLSWTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoD0starLLSWTest, uninitialized); \
    FRIEND_TEST(FFBtoD0starLLSWTest, clone);         \
    FRIEND_TEST(FFBtoD0starLLSWTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/LLSW/FFBtoD0starLLSW.hh"
#include "Hammer/Math/Constants.hh"
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

    TEST(FFBtoD0starLLSWTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FourMomentum pBmes{5279.63, 0., 0., 0.};
        double mD0star = 2349.0;
        double pz = 500.0;
        FourMomentum pD0starm{std::sqrt(mD0star * mD0star + pz * pz), 0., 0., pz};

        FFBtoD0starLLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, PID::BZERO), {Particle(pD0starm, PID::DSSD0STARMINUS)}, {});

        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BDSSD0STAR}, {0., -0.174548, 0.675519, 0.})};

        auto& t = ff.getTensor();

        EXPECT_NEAR(t.element({0}).real(), 0., 1e-10);
        EXPECT_NEAR(t.element({3}).real(), 0., 1e-10);

        for (IndexType idx = 1; idx <= 2; ++idx) {
            double comRe = compareVals(t.element({idx}).real(), ffEval.element({idx}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoD0starLLSWTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD0starLLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD0STARMINUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoD0starLLSWTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD0starLLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD0STARMINUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD0starLLSWTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starLLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5279.63, 0., 0., 0.};
        double mD0star = 2349.0;
        double pz = 500.0;
        FourMomentum pD0starm{std::sqrt(mD0star * mD0star + pz * pz), 0., 0., pz};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, PID::BZERO), {Particle(pD0starm, PID::DSSD0STARMINUS)}, {});

        ff.eval(Particle(pBmes, PID::BZERO), {Particle(pD0starm, PID::DSSD0STARMINUS)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFBtoD0starLLSWTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starLLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Leibovich:1997tu"));
        EXPECT_FALSE(set.checkReference("Leibovich:1997em"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Leibovich:1997tu"));
        EXPECT_TRUE(set.checkReference("Leibovich:1997em"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Leibovich:1997tu"));
        EXPECT_TRUE(set.checkReference("Leibovich:1997em"));

        const string tmpfile = testing::TempDir() + "hammer_llsw_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Leibovich:1997tu"), string::npos);
        EXPECT_NE(content.find("Leibovich:1997em"), string::npos);
    }

} // namespace Hammer
