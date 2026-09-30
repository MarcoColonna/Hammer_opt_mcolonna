///
/// @file  TestFFD2starDPiPW.cc
/// @brief Unit tests for the FFD2starDPiPW partial wave form factor class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFD2starDPiPWTest, evalDefaults);  \
    FRIEND_TEST(FFD2starDPiPWTest, evalCustomD);   \
    FRIEND_TEST(FFD2starDPiPWTest, uninitialized); \
    FRIEND_TEST(FFD2starDPiPWTest, clone)

#include "Hammer/FormFactors/PW/FFD2starDPiPW.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"

#include <complex>

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    // D2* (425, ~2461 MeV) -> D- (DMINUS=-411, ~1870 MeV) + pi+ (211, ~140 MeV)
    // D2* at rest; D and pi back-to-back along z.
    // p ~ 505.3 MeV/c, E_D- ~ 1936.7, E_pi+ ~ 524.3
    static const FourMomentum pD2star = []() noexcept { return FourMomentum{2461., 0., 0., 0.}; }();
    static const FourMomentum pDminus = []() noexcept { return FourMomentum{1936.7, 0., 0., 505.3}; }();
    static const FourMomentum pPi = []() noexcept { return FourMomentum{524.3, 0., 0., -505.3}; }();

    /// @brief Test that the default D=1+0i setting produces conj(D)=1 in the tensor.
    TEST(FFD2starDPiPWTest, evalDefaults) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD2starDPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pD2star, -PID::DSSD2STAR), {Particle(pDminus, PID::DMINUS), Particle(pPi, PID::PIPLUS)}, {});

        auto& t = ff.getTensor();

        // Default D=1+0i => conj(D)=1, stored at index 0
        EXPECT_DOUBLE_EQ(t.element({0}).real(), 1.0);
        EXPECT_DOUBLE_EQ(t.element({0}).imag(), 0.0);
    }

    /// @brief Test that a custom complex D value is conjugated and stored correctly.
    TEST(FFD2starDPiPWTest, evalCustomD) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD2starDPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        // Override D with a non-trivial complex value.
        // Settings path = prefix + group = "D**2*toDPi" + "PW" = "D**2*toDPiPW"
        set.changeSetting<complex<double>>("D**2*toDPiPW", "D", complex<double>(-0.2, 0.8));

        ff.calcUnits();
        ff.eval(Particle(pD2star, -PID::DSSD2STAR), {Particle(pDminus, PID::DMINUS), Particle(pPi, PID::PIPLUS)}, {});

        auto& t = ff.getTensor();

        // conj(-0.2 + 0.8i) = -0.2 - 0.8i
        EXPECT_NEAR(t.element({0}).real(), -0.2, 1e-12);
        EXPECT_NEAR(t.element({0}).imag(), -0.8, 1e-12);
    }

    TEST(FFD2starDPiPWTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFD2starDPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        ff.eval(Particle(pD2star, -PID::DSSD2STAR), {Particle(pDminus, PID::DMINUS), Particle(pPi, PID::PIPLUS)}, {});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFD2starDPiPWTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFD2starDPiPW ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pD2star, -PID::DSSD2STAR), {Particle(pDminus, PID::DMINUS), Particle(pPi, PID::PIPLUS)},
                     {});

        ff.eval(Particle(pD2star, -PID::DSSD2STAR), {Particle(pDminus, PID::DMINUS), Particle(pPi, PID::PIPLUS)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        EXPECT_NEAR(tOrig.element({0}).real(), tCloned.element({0}).real(), 1e-12);
        EXPECT_NEAR(tOrig.element({0}).imag(), tCloned.element({0}).imag(), 1e-12);
    }

} // namespace Hammer
