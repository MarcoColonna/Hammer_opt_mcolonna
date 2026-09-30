///
/// @file  TestFFBtoD1ISGW2.cc
/// @brief Tests for FFBtoD1ISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                    \
    FRIEND_TEST(FFBtoD1ISGW2Test, eval);          \
    FRIEND_TEST(FFBtoD1ISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoD1ISGW2Test, bsMode);        \
    FRIEND_TEST(FFBtoD1ISGW2Test, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoD1ISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoD1ISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD1mes{2.8240005676506059e+00, 3.4765542738957983e-01, 1.1098136867676383e+00,
                            8.7571280698419929e-01};
        FourMomentum kNuTau{1.0517172449734653e+00, -2.6820877777053831e-01, -5.5264138293147658e-01,
                            8.5367483081849993e-01};
        FourMomentum pTau{2.1983546059064443e+00, 1.4097211567917448e-01, -2.8907806513047413e-01,
                          1.2538877125849155e+00};

        // Evaluate the FF class
        FFBtoD1ISGW2 ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        set.changeSetting<bool>("BtoD**1ISGW2", "SmearQ2", false);
        // MeV momenta
        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD1mes, 10413), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        // Compare to direct EvtGen evaluation (also in GeV)
        double vv = -0.062548;
        double rr = -0.991994;
        double sppsm = -0.0721888;
        double spmsm = -0.0421309;
        double mb = pBmes.mass();  // GeV
        double mx = pD1mes.mass(); // GeV
        double sqmbmx = sqrt(mb * mx);
        double Fv1 = rr / sqmbmx;
        double Fv2 = mb * mb * sppsm / sqmbmx;
        double Fv3 = spmsm * sqmbmx;
        double Fa = 2 * vv * sqmbmx;
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSSD1}, {0, Fv1, Fv2, Fv3, Fa, 0, 0, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            double comRe = compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 2e-3); // Account for some EvtGen smearing
        }
    }

    TEST(FFBtoD1ISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1ISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5280.; // B0 (MeV)
        const double Mc = 2422.; // D1- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD1ISGW2Test, bsMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1ISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // Bs -> Ds1- signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5367.; // Bs (MeV)
        const double Mc = 2536.; // Ds1- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoD1ISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1ISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD1mes{2.8240005676506059e+00, 3.4765542738957983e-01, 1.1098136867676383e+00,
                            8.7571280698419929e-01};
        FourMomentum kNuTau{1.0517172449734653e+00, -2.6820877777053831e-01, -5.5264138293147658e-01,
                            8.5367483081849993e-01};
        FourMomentum pTau{2.1983546059064443e+00, 1.4097211567917448e-01, -2.8907806513047413e-01,
                          1.2538877125849155e+00};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(1000 * pBmes, -511),
                     {Particle(1000 * pD1mes, 10413), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD1mes, 10413), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
