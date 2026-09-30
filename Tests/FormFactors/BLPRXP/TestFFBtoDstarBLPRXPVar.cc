///
/// @file  TestFFBtoDstarBLPRXPVar.cc
/// @brief Tests for FFBtoDstarBLPRXPVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                           \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalBase);      \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalShifts);    \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalLOShifts);  \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalAsOM);      \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, uninitialized); \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, clone)

#define FFBLPRXPVARBASE_FRIENDS                          \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalBase);      \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalShifts);    \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalLOShifts);  \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, evalAsOM);      \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, uninitialized); \
    FRIEND_TEST(FFBtoDstarBLPRXPVarTest, clone)


#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLPRXP/FFBtoDstarBLPRXPVar.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    TEST(FFBtoDstarBLPRXPVarTest, evalBase) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoD*BLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoD*BLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXPVar", "WithHA1As2", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        array<double, 8> ffEval = {
            {-0.808529, 5393.47, 1.65619e-4, -7.0928e-06, -2.91102e-05, -3.24212e-08, 0.666864, -1.02329}};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            double comRe = compareVals(t.element({idx1, 0}).real(), ffEval[idx1]);
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDstarBLPRXPVarTest, evalShifts) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoD*BLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoD*BLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXPVar", "WithHA1As2", false);

        set.changeSetting<bool>("BtoD*BLPRXPVar", "WithAsOverM", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        double xi = ff.Xi((2 * 5280 * 2258.9451) / (2 * 5280 * 2007.));
        // Compare to direct evaluation
        array<array<double, 10>, 8> ffEval = {
            {{0., 0., -0.0517822, -0.00650038, 0.0124886, -0.206249, -0.0258911, 0.0140226, -0.0152338, 0.0127401},
             {0., 0., -220.492, -27.6791, -140.158, -63.8347, -8.01338, 12.9493, -13.9929, -142.98},
             {0., 0., -4.89456e-6, -6.1443e-7, -3.11127e-6, -0.0000197608, -2.48063e-6, -4.12895e-6, -3.1062e-7,
              -3.17392e-6},
             {0., 0., -0.000200533, -0.0000251735, -3.11127e-6, -0.0000433519, -5.44211e-6, 0.0000178189, -0.0000454497,
              -3.17392e-6},
             {0., 0., 0.0000927666, 0.0000116453, 3.11127e-6, 0.0000670935, 8.42246e-6, -0.0000102524, 0.0000205851,
              3.17392e-6},
             {0., 0., -2.68475e-8, -3.37026e-9, 0., -5.91469e-9, -7.4249e-10, -1.88504e-9, -6.19446e-9, 0.},
             {0., 0., -0.0160199, -0.00201103, -0.0101832, 0.14293, 0.0179425, -0.0312393, -0.00101666, -0.0103883},
             {0., 0., 0.0356667, 0.00447735, 0.0226718, -0.0657433, -0.00825296, 0.0173159, 0.00226349, 0.0231284}}};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 3; idx2 < 11; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), xi * ffEval[idx1][idx2 - 1]);
                EXPECT_NEAR(comRe, 1., 2e-3);
            }
        }
    }

    TEST(FFBtoDstarBLPRXPVarTest, evalLOShifts) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoD*BLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoD*BLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXPVar", "WithHA1As2", false);

        set.changeSetting<bool>("BtoD*BLPRXPVar", "WithAsOverM", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        double w = (2 * 5280 * 2258.9451) / (2 * 5280 * 2007.);

        // Compare to direct evaluation
        auto vLinXi = ff.XiLinearCoeffs(w);

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            double comRe = compareVals(t.element({idx1, 1}).real(), vLinXi[0] * t.element({idx1, 0}).real());
            EXPECT_NEAR(comRe, 1., 2e-3);
            comRe = compareVals(t.element({idx1, 2}).real(), vLinXi[1] * t.element({idx1, 0}).real());
            EXPECT_NEAR(comRe, 1., 2e-3);
        }
    }

    TEST(FFBtoDstarBLPRXPVarTest, evalAsOM) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDstarmes{2258.9451, -542.95362175877925474, -718.43449427445361377, 513.65086086722812801};
        FourMomentum kNuTau{763.45512706233784082, 59.818100429037400761, 266.26625745054465575, 713.01318784585764875};
        FourMomentum pTau{2257.5997729376621592, 483.13552132974185397, 452.16823682390895801, -1226.6640487130857768};

        // Evaluate the FF class
        FFBtoDstarBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoD*BLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoD*BLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoD*BLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoD*BLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoD*BLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoD*BLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoD*BLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        set.changeSetting<bool>("BtoD*BLPRXPVar", "WithHA1As2", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t = ff.getTensor();

        set.changeSetting<bool>("BtoD*BLPRXPVar", "WithAsOverM", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDstarmes, -413), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto t2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2258.9451) / (2 * 5280 * 2007.));
        // Compare to direct evaluation
        array<array<double, 10>, 8> ffEval = {
            {{0., 0., -0.0201334, -0.00252741, 0.0138751, -0.0116038, -0.00145666, 0., 0., 0.},
             {0., 0., -16.6377, -2.08858, -143.97, 2.51193, 0.315331, 0., 0., 0.},
             {0., 0., -9.07854e-7, -1.13966e-7, -3.53822e-6, 3.05842e-7, 3.83934e-8, 0., 0., 0.},
             {0., 0., -0.0000467594, -5.86986e-6, -3.41945e-6, 5.08693e-6, 6.38579e-7, 0., 0., 0.},
             {0., 0., 0.0000216847, 2.72215e-6, 3.18078e-6, -1.69822e-6, -2.13183e-7, 0., 0., 0.},
             {0., 0., -8.19453e-9, -1.02868e-9, 4.55397e-11, -3.02019e-11, -3.79134e-12, 0., 0., 0.},
             {0., 0., -0.00384439, -0.000482599, -0.0121355, -0.00674926, -0.000847256, 0., 0., 0.},
             {0., 0., 0.00555851, 0.000697777, 0.0251111, 0.00224833, 0.000282241, 0., 0., 0.}}};


        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 3; idx2 < 11; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real() - t2.element({idx1, idx2}).real(),
                                           xi * ffEval[idx1][idx2 - 1]);
                EXPECT_NEAR(comRe, 1., 2e-3);
            }
        }
    }

    TEST(FFBtoDstarBLPRXPVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBLPRXPVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSTAR);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_TRUE(std::isfinite(t.element({idx1, idx2}).real()))
                    << "element (" << idx1 << "," << idx2 << ") is not finite";
            }
        }
    }

    TEST(FFBtoDstarBLPRXPVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBLPRXPVar ff;
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

    TEST(FFBtoDstarBLPRXPVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarBLPRXPVar ff;
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
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

} // namespace Hammer
