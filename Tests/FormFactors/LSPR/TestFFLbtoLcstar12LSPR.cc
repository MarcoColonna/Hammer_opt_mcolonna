///
/// @file  TestFFLbtoLcstar12LSPR.cc
/// @brief Tests for FFLbtoLcstar12LSPR
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                          \
    FRIEND_TEST(FFLbtoLcstar12LSPRTest, eval);          \
    FRIEND_TEST(FFLbtoLcstar12LSPRTest, zeroRecoil);    \
    FRIEND_TEST(FFLbtoLcstar12LSPRTest, uninitialized); \
    FRIEND_TEST(FFLbtoLcstar12LSPRTest, clone);         \
    FRIEND_TEST(FFLbtoLcstar12LSPRTest, addRefs)

#include <filesystem>
#include <fstream>

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/LSPR/FFLbtoLcstar12LSPR.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFLbtoLcstar12LSPRTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3053.24065836299, 0, 0, -1608.80499684731};
        FourMomentum kNuTau{103.003318909325, -16.4797557212763, -28.5437742055744, -97.5876749994806};
        FourMomentum pTau{2463.75602272769, 16.4797557212763, 28.5437742055744, 1706.39267184679};

        // Evaluate the FF class
        FFLbtoLcstar12LSPR ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -14122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({12}, {FF_LBLCSTAR12},
                                               {0.800689720050897, 0.15180473106173, 0.164627941225841,
                                                -0.699710351366009, 0.00860587827688889, 0.811514188436004,
                                                -0.63695621301262, 0.00834708275883007, 0.854156249543776,
                                                -0.687821574518799, 0.00839946552333471, -0.0300044513272222})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFLbtoLcstar12LSPRTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcstar12LSPR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACSTAR12PLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 12; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFLbtoLcstar12LSPRTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcstar12LSPR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACSTAR12PLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFLbtoLcstar12LSPRTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcstar12LSPR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3053.24065836299, 0, 0, -1608.80499684731};
        FourMomentum kNuTau{103.003318909325, -16.4797557212763, -28.5437742055744, -97.5876749994806};
        FourMomentum pTau{2463.75602272769, 16.4797557212763, 28.5437742055744, 1706.39267184679};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pLbmes, -5122), {Particle(pLcmes, -14122), Particle(kNuTau, 16), Particle(pTau, -15)},
                     {});

        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -14122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 12; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFLbtoLcstar12LSPRTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcstar12LSPR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Papucci:2021pmj"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Papucci:2021pmj"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Papucci:2021pmj"));

        const string tmpfile = testing::TempDir() + "hammer_lspr_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Papucci:2021pmj"), string::npos);
    }

} // namespace Hammer
