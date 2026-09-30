///
/// @file  TestFFBtoD0starBLRVar.cc
/// @brief Tests for FFBtoD0starBLRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                         \
    FRIEND_TEST(FFBtoD0starBLRVarTest, eval);          \
    FRIEND_TEST(FFBtoD0starBLRVarTest, rename);        \
    FRIEND_TEST(FFBtoD0starBLRVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoD0starBLRVarTest, uninitialized); \
    FRIEND_TEST(FFBtoD0starBLRVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLR/FFBtoD0starBLRVar.hh"
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

    TEST(FFBtoD0starBLRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD0starmes{3117.27272727273, 0, 0, -2104.13622567512};
        FourMomentum kNuTau{49.97993502618, -24.8428481565870, -43.029075211928, -5.41530120089};
        FourMomentum pTau{2112.74733770109, 24.8428481565870, 43.029075211928, 2109.55152687601};

        // Evaluate the FF class
        FFBtoD0starBLRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pD0starmes, -10411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({4, 6}, {FF_BDSSD0STAR, FF_BDSSD0STAR_VAR},
                                               {0.495949833064184,
                                                0.708499761520263,
                                                0.164535722165732,
                                                -0.345774933396686,
                                                0.610104020438920,
                                                -0.478999981264232,
                                                -0.131663927760578,
                                                -0.188091325372254,
                                                -0.0436806668598272,
                                                0.259931083017855,
                                                0,
                                                0,
                                                0.720459193477596,
                                                1.02922741925371,
                                                0.239018678476773,
                                                0,
                                                1.71697794406059,
                                                -1.34801996952003,
                                                0.877081840967645,
                                                1.25297405852521,
                                                0.290979620278217,
                                                -0.146804930608091,
                                                1.71697794406059,
                                                -1.34801996952003})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-4);
            }
        }
    }


    TEST(FFBtoD0starBLRVarTest, rename) {

        SettingsHandler seth{};
        seth.addSetting<string>("Hammer", "Units", "MeV");
        seth.addSetting<vector<string>>("BtoD**0*BLRVar", "ErrNames", {"name1", "", "name2", "", ""});

        // Evaluate the FF class
        FFBtoD0starBLRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(seth);
        ff.initSettings();

        auto* names = seth.getNamedSettingValue<vector<string>>("BtoD**0*BLRVar", "ErrNames");
        EXPECT_EQ(names->size(), 5);
        EXPECT_EQ(names->at(0), "name1");
        EXPECT_EQ(names->at(1), "delta_ztp");
        EXPECT_EQ(names->at(2), "name2");
        EXPECT_EQ(names->at(3), "delta_chi1");
        EXPECT_EQ(names->at(4), "delta_chi2");
        auto* defaultNames = seth.getEntry("BtoD**0*BLRVar", "ErrNames", WTerm::COMMON)->getDefault<vector<string>>();
        EXPECT_EQ(defaultNames->size(), 5);
        EXPECT_EQ(defaultNames->at(0), "delta_zt1");
        EXPECT_EQ(defaultNames->at(1), "delta_ztp");
        EXPECT_EQ(defaultNames->at(2), "delta_zeta1");
        EXPECT_EQ(defaultNames->at(3), "delta_chi1");
        EXPECT_EQ(defaultNames->at(4), "delta_chi2");


        EXPECT_EQ(seth.getNamedSettingValue<double>("BtoD**0*BLRVar", "delta_zt1"), nullptr);
        EXPECT_EQ(*seth.getNamedSettingValue<double>("BtoD**0*BLRVar", "name1"), 0.0);
        EXPECT_EQ(*seth.getNamedSettingValue<double>("BtoD**0*BLRVar", "name2"), 0.0);
    }

    TEST(FFBtoD0starBLRVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD0starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD0STARPLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                EXPECT_TRUE(std::isfinite(t.element({idx1, idx2}).real()))
                    << "element (" << idx1 << "," << idx2 << ") is not finite";
            }
        }
    }

    TEST(FFBtoD0starBLRVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD0starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD0STARPLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD0starBLRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD0starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD0starmes{3117.27272727273, 0, 0, -2104.13622567512};
        FourMomentum kNuTau{49.97993502618, -24.8428481565870, -43.029075211928, -5.41530120089};
        FourMomentum pTau{2112.74733770109, 24.8428481565870, 43.029075211928, 2109.55152687601};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pD0starmes, -10411), Particle(kNuTau, 14), Particle(pTau, -13)},
                     {});

        ff.eval(Particle(pBmes, 511), {Particle(pD0starmes, -10411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 6; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }
} // namespace Hammer
