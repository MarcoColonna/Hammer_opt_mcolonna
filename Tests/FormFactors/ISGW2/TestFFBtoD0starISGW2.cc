///
/// @file  TestFFBtoD0starISGW2.cc
/// @brief Tests for FFBtoD0starISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                        \
    FRIEND_TEST(FFBtoD0starISGW2Test, eval);          \
    FRIEND_TEST(FFBtoD0starISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoD0starISGW2Test, bsMode);        \
    FRIEND_TEST(FFBtoD0starISGW2Test, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoD0starISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoD0starISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD0starmes{2.7206523538400131e+00, 3.5605736720613768e-01, 1.1540824843894137e+00,
                                7.9961099454269813e-01};
        FourMomentum kNuTau{1.1647033125190303e+00, -2.9539536290726093e-01, -6.1270578819281141e-01,
                            9.4544539920648796e-01};
        FourMomentum pTau{2.1887167521714703e+00, 1.5975676099933930e-01, -2.7328245749091501e-01,
                          1.2382189566384276e+00};

        // Evaluate the FF class
        FFBtoD0starISGW2 ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        // MeV momenta
        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD0starmes, 10411), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        // Compare to direct EvtGen evaluation (also in GeV)
        double uppum = -0.352846;
        double upmum = 0.759315;
        double mb = pBmes.mass();      // GeV
        double mx = pD0starmes.mass(); // GeV
        double gppgm = mb / sqrt(mb * mx) * uppum;
        double gpmgm = mx / sqrt(mb * mx) * upmum;
        double gp = (gppgm + gpmgm) / 2.0;
        double gm = (gppgm - gpmgm) / 2.0;
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BDSSD0STAR}, {0, gp, gm, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoD0starISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5280.; // B0 (MeV)
        const double Mc = 2349.; // D0*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD0starISGW2Test, bsMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // Bs -> Ds0*- signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5367.; // Bs (MeV)
        const double Mc = 2317.; // Ds0*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoD0starISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD0starmes{2.7206523538400131e+00, 3.5605736720613768e-01, 1.1540824843894137e+00,
                                7.9961099454269813e-01};
        FourMomentum kNuTau{1.1647033125190303e+00, -2.9539536290726093e-01, -6.1270578819281141e-01,
                            9.4544539920648796e-01};
        FourMomentum pTau{2.1887167521714703e+00, 1.5975676099933930e-01, -2.7328245749091501e-01,
                          1.2382189566384276e+00};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(1000 * pBmes, -511),
                     {Particle(1000 * pD0starmes, 10411), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD0starmes, 10411), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
