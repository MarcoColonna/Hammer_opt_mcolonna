///
/// @file  TestFFBtoDstarBGLX.cc
/// @brief Tests for FFBtoDstarBGLX
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFBtoDstarBGLXTest, eval);          \
    FRIEND_TEST(FFBtoDstarBGLXTest, evalOpt);          \
    FRIEND_TEST(FFBtoDstarBGLXTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDstarBGLXTest, uninitialized); \
    FRIEND_TEST(FFBtoDstarBGLXTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BGL/FFBtoDstarBGLX.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoDstarBGLXTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBGLX ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        //set.changeSetting<bool>("BtoD*BGLX", "ExtBranches", true);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSTAR},
                                               {-0.6845787512332282, 5057.250417739538, 0.00011113797765262371, 0.00022750456342981215,
                                                -0.00010514774368697381, 7.948449742654126e-11,0.4019876405768768,-0.9109755230439129})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; idx1++) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        };
    }

    TEST(FFBtoDstarBGLXTest, evalOpt) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBGLX ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<bool>("BtoD*BGLX", "OptZ", true);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSTAR},
                                               {-0.7921776782898936, 5616.443068618494,4.462384142511799e-6,0.00036546970467231854,
                                                -0.00014629821908523222,-4.234181740041228e-8,1.216137728315007,-1.0887944086920733})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; idx1++) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        };
    }
    
    TEST(FFBtoDstarBGLXTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBGLX ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSTARMINUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoDstarBGLXTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBGLX ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSTARMINUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDstarBGLXTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarBGLX ff;
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
