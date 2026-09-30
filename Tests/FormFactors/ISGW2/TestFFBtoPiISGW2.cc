///
/// @file  TestFFBtoPiISGW2.cc
/// @brief Tests for FFBtoPiISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                    \
    FRIEND_TEST(FFBtoPiISGW2Test, eval);          \
    FRIEND_TEST(FFBtoPiISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoPiISGW2Test, bsMode);        \
    FRIEND_TEST(FFBtoPiISGW2Test, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoPiISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoPiISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 = 19.8958 GeV^2)
        FourMomentum pBmes{5474.7423876404975e+00, 2844.8361598961946e-01, 1776.6650823948621e-01,
                           1409.3767837987914e+00};
        FourMomentum pPimes{7809.7394236330053e-01, -4349.7440268391885e-01, 5980.1941366927713e-01,
                            2088.3130366317992e-01};
        FourMomentum kNuMu{2952.7404598550526e+00, 2231.1045907528642e+00, -6669.6014486150135e-01,
                           3727.2426960535188e-01};
        FourMomentum pMu{1741.0279854221444e+00, -1511.6465720793257e+00, 2466.0723943171059e-01,
                         8278.2121053025981e-01};

        // Evaluate the FF class
        FFBtoPiISGW2 ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, -511), {Particle(pPimes, 211), Particle(kNuMu, -16), Particle(pMu, 15)}, {});

        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BPI}, {0., 0.417168, 1.22849, 0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 3; ++idx1) {
            double comRe = compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoPiISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5279.; // B+ (MeV)
        const double Mc = 135.;  // pi0 (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoPiISGW2Test, bsMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // Bs -> K- signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5367.; // Bs (MeV)
        const double Mc = 494.;  // K- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoPiISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5474.7423876404975e+00, 2844.8361598961946e-01, 1776.6650823948621e-01,
                           1409.3767837987914e+00};
        FourMomentum pPimes{7809.7394236330053e-01, -4349.7440268391885e-01, 5980.1941366927713e-01,
                            2088.3130366317992e-01};
        FourMomentum kNuMu{2952.7404598550526e+00, 2231.1045907528642e+00, -6669.6014486150135e-01,
                           3727.2426960535188e-01};
        FourMomentum pMu{1741.0279854221444e+00, -1511.6465720793257e+00, 2466.0723943171059e-01,
                         8278.2121053025981e-01};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, -511), {Particle(pPimes, 211), Particle(kNuMu, -16), Particle(pMu, 15)}, {});

        ff.eval(Particle(pBmes, -511), {Particle(pPimes, 211), Particle(kNuMu, -16), Particle(pMu, 15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
