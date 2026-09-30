///
/// @file  TestFFBtoDstarISGW2.cc
/// @brief Tests for FFBtoDstarISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                       \
    FRIEND_TEST(FFBtoDstarISGW2Test, eval);          \
    FRIEND_TEST(FFBtoDstarISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoDstarISGW2Test, bsMode);        \
    FRIEND_TEST(FFBtoDstarISGW2Test, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoDstarISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoDstarISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarISGW2 ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        //        //ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});
        //        TestableTensor FFres;
        //        dynamic_cast<Tensor&>(FFres) = ff.evalAtPSPoint({7.99026}, {5.280, 2.007});

        // Compare to direct evaluation
        // gammaji=-0.869422 mtx=2.15 r2=0.209807
        // The calculated form factors at t=7.99026 are ff=5.83968 gf=0.164454 af+=-0.133415 af-=0.152319
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        // Switched to eval on momenta (q^2 = 8.05199 GeV^2)
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSTAR},
                                               {0, 5841.5521663091869, 0.00016466092653151496, 0.00015250556021545483,
                                                -0.00013357097257008879, 0, 0, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 5; ++idx1) {
            double comRe = compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-3);
        }
    }

    TEST(FFBtoDstarISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5280.; // B0 (MeV)
        const double Mc = 2010.; // D*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDstarISGW2Test, bsMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // Bs -> Ds*- signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5367.; // Bs (MeV)
        const double Mc = 2112.; // Ds*- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoDstarISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
