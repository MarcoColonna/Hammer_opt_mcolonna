///
/// @file  TestFFBtoRhoOmegaISGW2.cc
/// @brief Tests for FFBtoRhoOmegaISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                          \
    FRIEND_TEST(FFBtoRhoOmegaISGW2Test, eval);          \
    FRIEND_TEST(FFBtoRhoOmegaISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoRhoOmegaISGW2Test, omegaMode);     \
    FRIEND_TEST(FFBtoRhoOmegaISGW2Test, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoRhoOmegaISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoRhoOmegaISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        // Momentum decs
        FourMomentum pBmes{5474.7423876404975e+00, 2844.8361598961946e-01, 1776.6650823948621e-01,
                           1409.3767837987914e+00};
        FourMomentum pRhomes{1613.2517093046810e+00, -7871.9937114550276e-01, 3761.4870333161349e-01,
                             1084.6451460327347e+00};
        FourMomentum kNuTau{1413.3906670454139e+00, 3690.9499116728806e-01, 1331.7693582553813e+00,
                            2963.6538529043410e-01};
        FourMomentum pTau{2448.1000112904030e+00, 7025.8799596783428e-01, -1530.2515533475087e+00,
                          2836.6252475622744e-02};

        // Evaluate the FF class
        FFBtoRhoOmegaISGW2 ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, -511), {Particle(pRhomes, 213), Particle(kNuTau, -16), Particle(pTau, 15)}, {});

        // Compare to direct evaluation
        // Ap,V,A0,A1,A12,T1,T2,T23
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BRHO}, {0, 1.15547, 0.487348, 0.384063, 0.22558, 0, 0, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 5; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoRhoOmegaISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoRhoOmegaISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5279.; // B+ (MeV)
        const double Mc = 775.;  // rho (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoRhoOmegaISGW2Test, omegaMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoRhoOmegaISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // B+ -> omega signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5279.; // B+ (MeV)
        const double Mc = 782.;  // omega (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoRhoOmegaISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoRhoOmegaISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5474.7423876404975e+00, 2844.8361598961946e-01, 1776.6650823948621e-01,
                           1409.3767837987914e+00};
        FourMomentum pRhomes{1613.2517093046810e+00, -7871.9937114550276e-01, 3761.4870333161349e-01,
                             1084.6451460327347e+00};
        FourMomentum kNuTau{1413.3906670454139e+00, 3690.9499116728806e-01, 1331.7693582553813e+00,
                            2963.6538529043410e-01};
        FourMomentum pTau{2448.1000112904030e+00, 7025.8799596783428e-01, -1530.2515533475087e+00,
                          2836.6252475622744e-02};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, -511), {Particle(pRhomes, 213), Particle(kNuTau, -16), Particle(pTau, 15)}, {});

        ff.eval(Particle(pBmes, -511), {Particle(pRhomes, 213), Particle(kNuTau, -16), Particle(pTau, 15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
