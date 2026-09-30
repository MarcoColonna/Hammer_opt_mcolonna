///
/// @file  TestFFD1DstarPiPW.cc
/// @brief Unit tests for the FFD1DstarPiPW partial wave form factor class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFD1DstarPiPWTest, evalDefaults);  \
    FRIEND_TEST(FFD1DstarPiPWTest, evalCustomSD);  \
    FRIEND_TEST(FFD1DstarPiPWTest, uninitialized); \
    FRIEND_TEST(FFD1DstarPiPWTest, clone)

#include "Hammer/FormFactors/PW/FFD1DstarPiPW.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"

#include <complex>

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    // D1 (10423, ~2421 MeV) -> D*- (DSTARMINUS=-413, ~2007 MeV) + pi+ (211, ~140 MeV)
    // D1 at rest; D* and pi back-to-back along z.
    // p ~ 356.3 MeV/c, E_D* ~ 2038.4, E_pi ~ 382.6
    static const FourMomentum pD1 = []() noexcept { return FourMomentum{2421., 0., 0., 0.}; }();
    static const FourMomentum pDstar = []() noexcept { return FourMomentum{2038.4, 0., 0., 356.3}; }();
    static const FourMomentum pPi = []() noexcept { return FourMomentum{382.6, 0., 0., -356.3}; }();

    /// @brief Test that the default S=0, D=1 settings produce the expected conjugated tensor values.
    TEST(FFD1DstarPiPWTest, evalDefaults) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD1DstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pD1, -PID::DSSD1), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)}, {});

        auto& t = ff.getTensor();

        // Default S=0+0i => conj(S)=0, stored at index 0
        EXPECT_DOUBLE_EQ(t.element({0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(t.element({0}).imag(), 0.0);

        // Default D=1+0i => conj(D)=1, stored at index 1
        EXPECT_DOUBLE_EQ(t.element({1}).real(), 1.0);
        EXPECT_DOUBLE_EQ(t.element({1}).imag(), 0.0);
    }

    /// @brief Test that custom complex S and D values are conjugated and stored correctly.
    TEST(FFD1DstarPiPWTest, evalCustomSD) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD1DstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        // Override S and D with non-trivial complex values.
        // Settings path = prefix + group = "D**1toD*Pi" + "PW" = "D**1toD*PiPW"
        set.changeSetting<complex<double>>("D**1toD*PiPW", "S", complex<double>(0.3, -0.4));
        set.changeSetting<complex<double>>("D**1toD*PiPW", "D", complex<double>(0.7, 0.2));

        ff.calcUnits();
        ff.eval(Particle(pD1, -PID::DSSD1), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)}, {});

        auto& t = ff.getTensor();

        // conj(0.3 - 0.4i) = 0.3 + 0.4i
        EXPECT_NEAR(t.element({0}).real(), 0.3, 1e-12);
        EXPECT_NEAR(t.element({0}).imag(), 0.4, 1e-12);

        // conj(0.7 + 0.2i) = 0.7 - 0.2i
        EXPECT_NEAR(t.element({1}).real(), 0.7, 1e-12);
        EXPECT_NEAR(t.element({1}).imag(), -0.2, 1e-12);
    }

    TEST(FFD1DstarPiPWTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFD1DstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        ff.eval(Particle(pD1, -PID::DSSD1), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)}, {});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFD1DstarPiPWTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD1DstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pD1, -PID::DSSD1), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)}, {});

        ff.eval(Particle(pD1, -PID::DSSD1), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 2; ++idx) {
            EXPECT_NEAR(tOrig.element({idx}).real(), tCloned.element({idx}).real(), 1e-12);
            EXPECT_NEAR(tOrig.element({idx}).imag(), tCloned.element({idx}).imag(), 1e-12);
        }
    }

} // namespace Hammer
