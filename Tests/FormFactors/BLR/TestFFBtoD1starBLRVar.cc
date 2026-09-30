///
/// @file  TestFFBtoD1starBLRVar.cc
/// @brief Tests for FFBtoD1starBLRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                         \
    FRIEND_TEST(FFBtoD1starBLRVarTest, eval);          \
    FRIEND_TEST(FFBtoD1starBLRVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoD1starBLRVarTest, uninitialized); \
    FRIEND_TEST(FFBtoD1starBLRVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLR/FFBtoD1starBLRVar.hh"
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

    TEST(FFBtoD1starBLRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD1starmes{3174.12206439394, 0, 0, -2045.65927751237};
        FourMomentum kNuTau{50.14709776263, -24.8428481565870, -43.029075211928, -6.78697251279};
        FourMomentum pTau{2055.73083784344, 24.8428481565870, 43.029075211928, 2052.44625002516};

        // Evaluate the FF class
        FFBtoD1starBLRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pD1starmes, -20413), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({8, 6}, {FF_BDSSD1STAR, FF_BDSSD1STAR_VAR},
                                               {0.719582925659466,
                                                1.02797560808495,
                                                0.208667593054135,
                                                0.223188464273391,
                                                -0.567249829058949,
                                                1.30912053298960,
                                                0.370603518630512,
                                                0.529433598043589,
                                                0.107469120586971,
                                                0.112086113148827,
                                                -0.174620874871706,
                                                0,
                                                0.365643880996805,
                                                0.522348401424007,
                                                0.106030904628056,
                                                0.567249829058949,
                                                0,
                                                -0.567249829058949,
                                                -1.18565830953627,
                                                -1.69379758505181,
                                                -0.343822034699932,
                                                -0.615817418460634,
                                                0.567249829058949,
                                                -0.567249829058949,
                                                0.782229873159499,
                                                1.11747124737071,
                                                0.226834210522220,
                                                0.0485675894016848,
                                                -0.567249829058949,
                                                0,
                                                -0.668353566591224,
                                                -0.954790809416035,
                                                -0.193811894469182,
                                                -0.0485675894016848,
                                                0.567249829058949,
                                                0,
                                                0.862587515435831,
                                                1.23226787919404,
                                                0.250136647530353,
                                                -0.0485675894016848,
                                                -0.567249829058949,
                                                0,
                                                0.392418134101618,
                                                0.560597334430882,
                                                0.113795011796223,
                                                0.567249829058949,
                                                0,
                                                0.567249829058949})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-4);
            }
        }
    }

    TEST(FFBtoD1starBLRVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD1starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD1STARPLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                EXPECT_TRUE(std::isfinite(t.element({idx1, idx2}).real()))
                    << "element (" << idx1 << "," << idx2 << ") is not finite";
            }
        }
    }

    TEST(FFBtoD1starBLRVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD1starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD1STARPLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD1starBLRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD1starmes{3174.12206439394, 0, 0, -2045.65927751237};
        FourMomentum kNuTau{50.14709776263, -24.8428481565870, -43.029075211928, -6.78697251279};
        FourMomentum pTau{2055.73083784344, 24.8428481565870, 43.029075211928, 2052.44625002516};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pD1starmes, -20413), Particle(kNuTau, 14), Particle(pTau, -13)},
                     {});

        ff.eval(Particle(pBmes, 511), {Particle(pD1starmes, -20413), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

} // namespace Hammer
