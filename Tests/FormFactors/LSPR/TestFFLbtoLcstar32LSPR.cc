///
/// @file  TestFFLbtoLcstar32LSPR.cc
/// @brief Tests for FFLbtoLcstar32LSPR
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                          \
    FRIEND_TEST(FFLbtoLcstar32LSPRTest, eval);          \
    FRIEND_TEST(FFLbtoLcstar32LSPRTest, zeroRecoil);    \
    FRIEND_TEST(FFLbtoLcstar32LSPRTest, uninitialized); \
    FRIEND_TEST(FFLbtoLcstar32LSPRTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/LSPR/FFLbtoLcstar32LSPR.hh"
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

    TEST(FFLbtoLcstar32LSPRTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3067.17304270463, 0, 0, -1586.48210638947};
        FourMomentum kNuTau{103.858273281499, -16.4797557212763, -28.5437742055744, -98.4896519157352};
        FourMomentum pTau{2448.96868401388, 16.4797557212763, 28.5437742055744, 1684.97175830520};

        // Evaluate the FF class
        FFLbtoLcstar32LSPR ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4124), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval",
                      MD::makeVector({16}, {FF_LBLCSTAR32},
                                     {0.720401386256984, 1.24417777068813, 1.26182697446594, -0.470732504336548,
                                      -0.0170396347149046, 0.0207258015212351, 0.738050590034792, -0.511494766669524,
                                      0.0281547899021519, -0.0207258015212351, 0.781633411513197, -0.496146403812827,
                                      0.0276299806183191, -0.0207258015212351, 0, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 16; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFLbtoLcstar32LSPRTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcstar32LSPR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACSTAR32PLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 16; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFLbtoLcstar32LSPRTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcstar32LSPR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACSTAR32PLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFLbtoLcstar32LSPRTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcstar32LSPR ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{3067.17304270463, 0, 0, -1586.48210638947};
        FourMomentum kNuTau{103.858273281499, -16.4797557212763, -28.5437742055744, -98.4896519157352};
        FourMomentum pTau{2448.96868401388, 16.4797557212763, 28.5437742055744, 1684.97175830520};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4124), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4124), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 16; ++idx) {
            EXPECT_NEAR(tOrig.element({idx}).real(), tCloned.element({idx}).real(), 1e-12);
        }
    }


} // namespace Hammer
