///
/// @file  TestFFBtoDstarCLN.cc
/// @brief Tests for FFBtoDstarCLN
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFBtoDstarCLNTest, eval);          \
    FRIEND_TEST(FFBtoDstarCLNTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDstarCLNTest, uninitialized); \
    FRIEND_TEST(FFBtoDstarCLNTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/CLN/FFBtoDstarCLN.hh"
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

    TEST(FFBtoDstarCLNTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarCLN ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSTAR},
                                               {-0.6834702353581162, 5426.573911221199, 0.00016704618345474333,
                                                0.0001326215559184461, -0.00010442321895879478, 0., 0., 0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 5; ++idx1) {
            // std::cout<<t.element({idx1}).real()<<'\n';
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoDstarCLNTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarCLN ff;
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

    TEST(FFBtoDstarCLNTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarCLN ff;
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

    TEST(FFBtoDstarCLNTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarCLN ff;
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
