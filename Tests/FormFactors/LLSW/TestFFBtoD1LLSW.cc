///
/// @file  TestFFBtoD1LLSW.cc
/// @brief Unit tests for the FFBtoD1LLSW form factor class (B -> D1)
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                   \
    FRIEND_TEST(FFBtoD1LLSWTest, eval);          \
    FRIEND_TEST(FFBtoD1LLSWTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoD1LLSWTest, uninitialized); \
    FRIEND_TEST(FFBtoD1LLSWTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/LLSW/FFBtoD1LLSW.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoD1LLSWTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FourMomentum pBmes{5279.63, 0., 0., 0.};
        double mD1 = 2423.0;
        double pz = 500.0;
        FourMomentum pD1m{std::sqrt(mD1 * mD1 + pz * pz), 0., 0., pz};

        FFBtoD1LLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, PID::BZERO), {Particle(pD1m, PID::DSSD1MINUS)}, {});

        Tensor ffEval{"ffEval",
                      MD::makeVector({8}, {FF_BDSSD1}, {0., -0.386631, -2.08835, 1.17527, -0.743488, 0., 0., 0.})};

        auto& t = ff.getTensor();

        EXPECT_NEAR(t.element({0}).real(), 0., 1e-10);
        for (IndexType idx = 5; idx <= 7; ++idx) {
            EXPECT_NEAR(t.element({idx}).real(), 0., 1e-10);
        }

        for (IndexType idx = 1; idx <= 4; ++idx) {
            double comRe = compareVals(t.element({idx}).real(), ffEval.element({idx}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoD1LLSWTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD1LLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD1MINUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoD1LLSWTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD1LLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD1MINUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD1LLSWTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1LLSW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5279.63, 0., 0., 0.};
        double mD1 = 2423.0;
        double pz = 500.0;
        FourMomentum pD1m{std::sqrt(mD1 * mD1 + pz * pz), 0., 0., pz};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, PID::BZERO), {Particle(pD1m, PID::DSSD1MINUS)}, {});

        ff.eval(Particle(pBmes, PID::BZERO), {Particle(pD1m, PID::DSSD1MINUS)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
