///
/// @file  TestFFLbtoLcBLRSVar.cc
/// @brief Tests for FFLbtoLcBLRSVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                       \
    FRIEND_TEST(FFLbtoLcBLRSVarTest, eval);          \
    FRIEND_TEST(FFLbtoLcBLRSVarTest, zeroRecoil);    \
    FRIEND_TEST(FFLbtoLcBLRSVarTest, uninitialized); \
    FRIEND_TEST(FFLbtoLcBLRSVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLRS/FFLbtoLcBLRSVar.hh"
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

    TEST(FFLbtoLcBLRSVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{2514.6000000000004, 380.2511160960518, -878.708483534587, -425.0854616098454};
        FourMomentum pTau{2456.13415321897, -581.6942711749232, 1418.910848905593, 723.6624190525239};
        FourMomentum kNuTau{649.2658467810297, 201.44315507887137, -540.2023653710064, -298.5769574426786};

        // Evaluate the FF class
        FFLbtoLcBLRSVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({12, 5}, {FF_LBLC, FF_LBLC_VAR},
                                               {0.7312555660768754,
                                                0.09005025183411991,
                                                0.0045025125917059995,
                                                1.1624918768285328e-7,
                                                0.,
                                                1.1660722177550018,
                                                0.1435956206514121,
                                                0.007179781032570612,
                                                1.1624918768285328e-7,
                                                -1.1624918768285328e-7,
                                                1.195461806802338,
                                                0.14721479295968237,
                                                0.0073607396479841245,
                                                1.1624918768285328e-7,
                                                -1.1624918768285328e-7,
                                                -0.31472386413207737,
                                                -0.038756577779431874,
                                                -0.0019378288889715953,
                                                0.,
                                                1.1624918768285328e-7,
                                                -0.0851912009508891,
                                                -0.01049084541104497,
                                                -0.0005245422705522489,
                                                0.,
                                                0.,
                                                0.7527695959975337,
                                                0.09269959072765581,
                                                0.004634979536382794,
                                                1.1624918768285328e-7,
                                                0.,
                                                -0.3856860143396382,
                                                -0.04749519091733134,
                                                -0.002374759545866569,
                                                0.,
                                                1.1624918768285328e-7,
                                                0.11604815859155679,
                                                0.014290716393613024,
                                                0.0007145358196806519,
                                                0.,
                                                0.,
                                                0.8056467077825589,
                                                0.09921112712258268,
                                                0.004960556356129138,
                                                1.1624918768285328e-7,
                                                0.,
                                                -0.37394516263654165,
                                                -0.04604936718393667,
                                                -0.0023024683591968357,
                                                0.,
                                                1.1624918768285328e-7,
                                                0.11780903740914717,
                                                0.014507559298240859,
                                                0.0007253779649120436,
                                                0.,
                                                0.,
                                                0.013524178395048072,
                                                0.0016654309749152896,
                                                0.00008327154874576455,
                                                0.,
                                                0.})};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFLbtoLcBLRSVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcBLRSVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACPLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                EXPECT_TRUE(std::isfinite(t.element({idx1, idx2}).real()))
                    << "element (" << idx1 << "," << idx2 << ") is not finite";
            }
        }
    }

    TEST(FFLbtoLcBLRSVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFLbtoLcBLRSVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::LAMBDAB);
        const double Mc = pdg.getMass(PID::LAMBDACPLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFLbtoLcBLRSVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFLbtoLcBLRSVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pLbmes{5620., 0., 0., 0.};
        FourMomentum pLcmes{2514.6000000000004, 380.2511160960518, -878.708483534587, -425.0854616098454};
        FourMomentum kNuTau{649.2658467810297, 201.44315507887137, -540.2023653710064, -298.5769574426786};
        FourMomentum pTau{2456.13415321897, -581.6942711749232, 1418.910848905593, 723.6624190525239};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        ff.eval(Particle(pLbmes, -5122), {Particle(pLcmes, -4122), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 12; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }


} // namespace Hammer
