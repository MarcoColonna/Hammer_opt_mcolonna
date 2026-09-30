///
/// @file  TestFFD1starDstarPiPW.cc
/// @brief Unit tests for the FFD1starDstarPiPW partial wave form factor class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                         \
    FRIEND_TEST(FFD1starDstarPiPWTest, evalDefaults);  \
    FRIEND_TEST(FFD1starDstarPiPWTest, evalCustomS);   \
    FRIEND_TEST(FFD1starDstarPiPWTest, uninitialized); \
    FRIEND_TEST(FFD1starDstarPiPWTest, clone)

#include "Hammer/FormFactors/PW/FFD1starDstarPiPW.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"

#include <complex>

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    // D1* (20423, ~2427 MeV) -> D*- (DSTARMINUS=-413, ~2007 MeV) + pi+ (211, ~140 MeV)
    // D1* at rest; D* and pi back-to-back along z.
    // p ~ 361.7 MeV/c, E_D* ~ 2039.3, E_pi ~ 387.7
    static const FourMomentum pD1star = []() noexcept { return FourMomentum{2427., 0., 0., 0.}; }();
    static const FourMomentum pDstar = []() noexcept { return FourMomentum{2039.3, 0., 0., 361.7}; }();
    static const FourMomentum pPi = []() noexcept { return FourMomentum{387.7, 0., 0., -361.7}; }();

    /// @brief Test that the default S=1+0i setting produces conj(S)=1 in the tensor.
    TEST(FFD1starDstarPiPWTest, evalDefaults) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD1starDstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pD1star, -PID::DSSD1STAR), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)},
                {});

        auto& t = ff.getTensor();

        // Default S=1+0i => conj(S)=1, stored at index 0
        EXPECT_DOUBLE_EQ(t.element({0}).real(), 1.0);
        EXPECT_DOUBLE_EQ(t.element({0}).imag(), 0.0);
    }

    /// @brief Test that a custom complex S value is conjugated and stored correctly.
    TEST(FFD1starDstarPiPWTest, evalCustomS) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD1starDstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        // Override S with a non-trivial complex value.
        // Settings path = prefix + group = "D**1*toD*Pi" + "PW" = "D**1*toD*PiPW"
        set.changeSetting<complex<double>>("D**1*toD*PiPW", "S", complex<double>(0.5, -0.6));

        ff.calcUnits();
        ff.eval(Particle(pD1star, -PID::DSSD1STAR), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)},
                {});

        auto& t = ff.getTensor();

        // conj(0.5 - 0.6i) = 0.5 + 0.6i
        EXPECT_NEAR(t.element({0}).real(), 0.5, 1e-12);
        EXPECT_NEAR(t.element({0}).imag(), 0.6, 1e-12);
    }

    TEST(FFD1starDstarPiPWTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFD1starDstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        ff.eval(Particle(pD1star, -PID::DSSD1STAR), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)},
                {});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFD1starDstarPiPWTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD1starDstarPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pD1star, -PID::DSSD1STAR),
                     {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)}, {});

        ff.eval(Particle(pD1star, -PID::DSSD1STAR), {Particle(pDstar, PID::DSTARMINUS), Particle(pPi, PID::PIPLUS)},
                {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        EXPECT_NEAR(tOrig.element({0}).real(), tCloned.element({0}).real(), 1e-12);
        EXPECT_NEAR(tOrig.element({0}).imag(), tCloned.element({0}).imag(), 1e-12);
    }

} // namespace Hammer
