///
/// @file  TestFFBtoD2starISGW2.cc
/// @brief Tests for FFBtoD2starISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                        \
    FRIEND_TEST(FFBtoD2starISGW2Test, eval);          \
    FRIEND_TEST(FFBtoD2starISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoD2starISGW2Test, bsMode);        \
    FRIEND_TEST(FFBtoD2starISGW2Test, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoD2starISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoD2starISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD2starmes{2.8848738200870283e+00, -5.2810511065387389e-01, 1.0249552509378312e+00,
                                9.7358396705915939e-01};
        FourMomentum kNuTau{6.3239441837503607e-01, 2.3275878019939589e-01, 2.2488816988188981e-01,
                            5.4329675287006918e-01};
        FourMomentum pTau{2.5568041800684504e+00, 5.1576509575269402e-01, -9.8174918211403384e-01,
                          1.4663946304583855e+00};

        // Evaluate the FF class
        FFBtoD2starISGW2 ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        // MeV momenta
        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD2starmes, 415), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        // Compare to direct EvtGen evaluation (also in GeV)
        double hf = 0.0130539;
        double kf = 0.647262;
        double bppbm = 0.000378403;
        double bpmbm = -0.0237689;
        double mb = pBmes.mass();      // GeV
        double mx = pD2starmes.mass(); // GeV
        double sqmbmx = sqrt(mb * mx);
        double Ka1 = kf * mb / sqmbmx;
        double Ka2 = mb * mb * mb * bppbm / sqmbmx;
        double Ka3 = bpmbm * mb * sqmbmx;
        double Kv = 2 * hf * mb * sqmbmx;
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSSD2STAR}, {0, Ka1, Ka2, Ka3, Kv, 0, 0, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoD2starISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD2starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5280.; // B0 (MeV)
        const double Mc = 2460.; // D2*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD2starISGW2Test, bsMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD2starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // Bs -> Ds2*- signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5367.; // Bs (MeV)
        const double Mc = 2573.; // Ds2*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoD2starISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD2starISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{6.0740724185305135e+00, 2.2041876529821600e-01, 2.6809423870568733e-01,
                           2.9832753503876139e+00};
        FourMomentum pD2starmes{2.8848738200870283e+00, -5.2810511065387389e-01, 1.0249552509378312e+00,
                                9.7358396705915939e-01};
        FourMomentum kNuTau{6.3239441837503607e-01, 2.3275878019939589e-01, 2.2488816988188981e-01,
                            5.4329675287006918e-01};
        FourMomentum pTau{2.5568041800684504e+00, 5.1576509575269402e-01, -9.8174918211403384e-01,
                          1.4663946304583855e+00};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(1000 * pBmes, -511),
                     {Particle(1000 * pD2starmes, 415), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        ff.eval(Particle(1000 * pBmes, -511),
                {Particle(1000 * pD2starmes, 415), Particle(1000 * kNuTau, -16), Particle(1000 * pTau, 15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
