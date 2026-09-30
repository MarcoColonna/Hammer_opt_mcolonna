///
/// @file  TestFFBtoD1BLRVar.cc
/// @brief Tests for FFBtoD1BLRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFBtoD1BLRVarTest, eval);          \
    FRIEND_TEST(FFBtoD1BLRVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoD1BLRVarTest, uninitialized); \
    FRIEND_TEST(FFBtoD1BLRVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLR/FFBtoD1BLRVar.hh"
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

    TEST(FFBtoD1BLRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD1mes{3171.36751893939, 0, 0, -2048.49479867141};
        FourMomentum kNuTau{50.13801061207, -24.8428481565870, -43.029075211928, -6.71950072636};
        FourMomentum pTau{2058.49447044853, 24.8428481565870, 43.029075211928, 2055.21429939777};

        // Evaluate the FF class
        FFBtoD1BLRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pD1mes, -10413), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({8, 8}, {FF_BDSSD1, FF_BDSSD1_VAR},
                                               {-0.739944645300832,
                                                -1.05706377900119,
                                                -0.454953279247777,
                                                0.157692870471125,
                                                0.0435629828789301,
                                                -1.52410523149317,
                                                -0.157460974963249,
                                                0.254017538582194,
                                                -0.843405995077826,
                                                -1.20486570725404,
                                                -0.518566254428249,
                                                0.157728807451071,
                                                -0.107147011625735,
                                                0.0787304874816245,
                                                0,
                                                0.118095731222437,
                                                -1.13030094119333,
                                                -1.61471563027618,
                                                -0.694962958376009,
                                                -0.0670443728172383,
                                                -0.229039764511052,
                                                -0.549835520309406,
                                                -0.0681666609168640,
                                                0.274917760154703,
                                                0.905362519914457,
                                                1.29337502844922,
                                                0.556660082560165,
                                                -0.00123823585356361,
                                                0.275424933988451,
                                                -0.803853058891600,
                                                -0.0681666609168640,
                                                -0.106108547718589,
                                                -0.643238103665703,
                                                -0.918911576665291,
                                                -0.395493482462943,
                                                0.0682826086708020,
                                                -0.0463851694773990,
                                                0.254017538582194,
                                                0,
                                                0.381026307873291,
                                                0.0249154323981526,
                                                0.0355934748545037,
                                                0.0153191968418223,
                                                0.0339673827044941,
                                                -0.0558648218978970,
                                                -0.254017538582194,
                                                0,
                                                -0.381026307873291,
                                                -0.660647522031832,
                                                -0.943782174331188,
                                                -0.406197623679137,
                                                0.0339673827044941,
                                                -0.0558648218978970,
                                                0.254017538582194,
                                                0,
                                                0.381026307873291,
                                                -0.203766349305352,
                                                -0.291094784721932,
                                                -0.125285275602107,
                                                -0.399190261202133,
                                                -0.320795755798354,
                                                0.549835520309406,
                                                0.0681666609168640,
                                                -0.274917760154703})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-4);
            }
        }
    }

    TEST(FFBtoD1BLRVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD1BLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD1PLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                EXPECT_TRUE(std::isfinite(t.element({idx1, idx2}).real()))
                    << "element (" << idx1 << "," << idx2 << ") is not finite";
            }
        }
    }

    TEST(FFBtoD1BLRVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD1BLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD1PLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD1BLRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD1BLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD1mes{3171.36751893939, 0, 0, -2048.49479867141};
        FourMomentum kNuTau{50.13801061207, -24.8428481565870, -43.029075211928, -6.71950072636};
        FourMomentum pTau{2058.49447044853, 24.8428481565870, 43.029075211928, 2055.21429939777};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pD1mes, -10413), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pD1mes, -10413), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

} // namespace Hammer
