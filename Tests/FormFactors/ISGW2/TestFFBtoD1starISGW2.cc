///
/// @file  TestFFBtoD1starISGW2.cc
/// @brief Tests for FFBtoD1starISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                        \
    FRIEND_TEST(FFBtoD1starISGW2Test, eval);          \
    FRIEND_TEST(FFBtoD1starISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoD1starISGW2Test, bsMode);        \
    FRIEND_TEST(FFBtoD1starISGW2Test, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoD1starISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoD1starISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD1starmes{2.8324053852172666e+00, 3.4694469299009156e-01, 1.1061032382391849e+00,
                                8.8195903786572871e-01};
        FourMomentum kNuTau{1.0425915030586572e+00, -2.6601336043581686e-01, -5.4778975227448157e-01,
                            8.4626256068835759e-01};
        FourMomentum pTau{2.1990755302545901e+00, 1.3948743274394138e-01, -2.9021924725901632e-01,
                          1.2550537518335279e+00};

        // Evaluate the FF class
        FFBtoD1starISGW2 ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        set.changeSetting<bool>("BtoD**1*ISGW2", "SmearQ2", false);
        // MeV momenta
        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD1starmes, 20413), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        // Compare to direct EvtGen evaluation (also in GeV)
        double ql = 0.0665166;
        double ll = 0.205976;
        double cppcm = 0.00451374;
        double cpmcm = -0.154666;
        double mb = pBmes.mass();      // GeV
        double mx = pD1starmes.mass(); // GeV
        double sqmbmx = sqrt(mb * mx);
        double Gv1 = ll / sqmbmx;
        double Gv2 = mb * mb * cppcm / sqmbmx;
        double Gv3 = cpmcm * sqmbmx;
        double Ga = 2 * ql * sqmbmx;
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSSD1STAR}, {0, Gv1, Gv2, Gv3, Ga, 0, 0, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoD1starISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5280.; // B0 (MeV)
        const double Mc = 2420.; // D1*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD1starISGW2Test, bsMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // Bs -> Ds1*- signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5367.; // Bs (MeV)
        const double Mc = 2460.; // Ds1*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoD1starISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD1starmes{2.8324053852172666e+00, 3.4694469299009156e-01, 1.1061032382391849e+00,
                                8.8195903786572871e-01};
        FourMomentum kNuTau{1.0425915030586572e+00, -2.6601336043581686e-01, -5.4778975227448157e-01,
                            8.4626256068835759e-01};
        FourMomentum pTau{2.1990755302545901e+00, 1.3948743274394138e-01, -2.9021924725901632e-01,
                          1.2550537518335279e+00};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(1000 * pBmes, -511),
                     {Particle(1000 * pD1starmes, 20413), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD1starmes, 20413), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
