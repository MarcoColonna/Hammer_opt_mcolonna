///
/// @file  TestFFBtoDstarBLPRXP.cc
/// @brief Tests for FFBtoDstarBLPRXP
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                        \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, AsOverM);       \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, OneOverMbMc);   \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, OneOverMb2);    \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, eval);          \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, evalAs2);       \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, uninitialized); \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, clone)

#define FFBLPRXPBASE_FRIENDS                          \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, AsOverM);       \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, OneOverMbMc);   \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, OneOverMb2);    \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, eval);          \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, evalAs2);       \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, uninitialized); \
    FRIEND_TEST(FFBtoDstarBLPRXPTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLPRXP/FFBLPRXPBase.hh"
#include "Hammer/FormFactors/BLPRXP/FFBtoDstarBLPRXP.hh"
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

    TEST(FFBtoDstarBLPRXPTest, AsOverM) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXP", "la2", 0.13);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXP", "WithHA1As2", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});
        auto tDiff = ff.getTensor();

        set.changeSetting<bool>("BtoD*BLPRXP", "WithAsOverM", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto tDiff2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2258.9451) / (2 * 5280 * 2007.));
        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{
            "ffEval",
            MD::makeVector({8}, {FF_BDSTAR},
                           {xi * (-0.616534) * 0.0394438, xi * 6919.24 * (-0.00479026), xi * 0.000153596 * 0.0147054,
                            xi * ((-0.0000583839) * 0.112646 + 0.000153596 * (-0.137499)),
                            xi * ((-0.0000583839) * 0.112646 + (-0.000153596) * (-0.137499)),
                            xi * 2.90901e-8 * (-0.247977), xi * (0.502719 * (-0.00874116) + (-1.11925) * (-0.0151578)),
                            xi * ((-1.11925) * (-0.00874116) + 0.502719 * (-0.0151578))})};

        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            double comRe = compareVals(tDiff.element({idx1}).real() - tDiff2.element({idx1}).real(),
                                       ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDstarBLPRXPTest, OneOverMbMc) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXP", "la2", 0.13);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXP", "WithHA1As2", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});
        auto tDiff = ff.getTensor();

        set.changeSetting<bool>("BtoD*BLPRXP", "With1OverMbMc", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto tDiff2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2258.9451) / (2 * 5280 * 2007.));
        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{
            "ffEval",
            MD::makeVector({8}, {FF_BDSTAR},
                           {xi * (-0.616534) * (-0.000391075), xi * 6919.24 * 0.00237044, xi * 0.000153596 * 0.00237044,
                            xi * ((-0.0000583839) * (-0.0157677) + 0.000153596 * 0.0181381),
                            xi * ((-0.0000583839) * (-0.0157677) + (-0.000153596) * 0.0181381),
                            xi * 2.90901e-8 * (-0.0157677), xi * (0.502719 * (-0.00237044) + (-1.11925) * 0.),
                            xi * ((-1.11925) * (-0.00237044) + 0.502719 * 0.)})};

        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            double comRe = compareVals(tDiff.element({idx1}).real() - tDiff2.element({idx1}).real(),
                                       ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDstarBLPRXPTest, OneOverMb2) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXP", "la2", 0.13);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXP", "WithHA1As2", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});
        auto tDiff = ff.getTensor();

        set.changeSetting<bool>("BtoD*BLPRXP", "With1OverMb2", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto tDiff2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2258.9451) / (2 * 5280 * 2007.));
        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval",
                      MD::makeVector({8}, {FF_BDSTAR},
                                     {xi * (-0.616534) * (-0.000558246), xi * 6919.24 * 0.00059352,
                                      xi * 0.000153596 * (-0.000558246),
                                      xi * ((-0.0000583839) * 0. + 0.000153596 * (-0.000558246)),
                                      xi * ((-0.0000583839) * 0. + (-0.000153596) * (-0.000558246)),
                                      xi * 2.90901e-8 * 0., xi * (0.502719 * 0.000665813 + (-1.11925) * (-0.00122406)),
                                      xi * ((-1.11925) * 0.000665813 + 0.502719 * (-0.00122406))})};

        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            double comRe = compareVals(tDiff.element({idx1}).real() - tDiff2.element({idx1}).real(),
                                       ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDstarBLPRXPTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXP", "la2", 0.13);
        set.changeSetting<double>("BtoD*BLPRXP", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoD*BLPRXP", "cSt", 1.55);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXP", "WithHA1As2", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BDSTAR},
                                               {-0.808529, 5393.47, 1.65619e-4, -7.0928e-06, -2.91102e-05, -3.24212e-08,
                                                0.666864, -1.02329})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoDstarBLPRXPTest, evalAs2) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXP", "la2", 0.13);
        set.changeSetting<double>("BtoD*BLPRXP", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoD*BLPRXP", "cSt", 1.55);
        ff.calcUnits();

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t = ff.getTensor();

        set.changeSetting<bool>("BtoD*BLPRXP", "WithHA1As2", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2258.9451) / (2 * 5280 * 2007.));

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval",
                      MD::makeVector({8}, {FF_BDSTAR}, {0., xi * 6919.24 * (-0.00929691), 0., 0., 0., 0., 0., 0.})};

        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            double comRe =
                compareVals(t.element({idx1}).real() - t2.element({idx1}).real(), ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDstarBLPRXPTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBLPRXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSTAR);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoDstarBLPRXPTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBLPRXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSTAR);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1.;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDstarBLPRXPTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarBLPRXP ff;
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
