///
/// @file  TestFFBtoDBLPRXPVar.cc
/// @brief Tests for FFBtoDBLPRXPVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                       \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalBase);      \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalShift);     \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalLOShift);   \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalAsOM);      \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, uninitialized); \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, clone)

#define FFBLPRXPVARBASE_FRIENDS                      \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalBase);      \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalShift);     \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalLOShift);   \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, evalAsOM);      \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, uninitialized); \
    FRIEND_TEST(FFBtoDBLPRXPVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLPRXP/FFBtoDBLPRXPVar.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    TEST(FFBtoDBLPRXPVarTest, evalBase) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoDBLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoDBLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        array<double, 4> ffEval = {{5228.31, 0.783358, 0.8981, 0.0001471}};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            double comRe = compareVals(t.element({idx1, 0}).real(), ffEval[idx1]);
            EXPECT_NEAR(comRe, 1., 2e-3);
        }
    }

    TEST(FFBtoDBLPRXPVarTest, evalShift) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoDBLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoDBLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        set.changeSetting<bool>("BtoDBLPRXPVar", "WithAsOverM", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        double xi = ff.Xi((2 * 5280 * 2381.3371929065643) / (2. * 5280. * 1869));
        // Compare to direct evaluation
        array<array<double, 10>, 4> ffEval = {
            {{0., 0., -2304.46, -631.708, 6913.39, -557.107, -152.716, -77.2002, -448.553, 1345.66},
             {0., 0., -0.322348, -0.0883631, 0.967043, -0.0773993, -0.021217, -0.0132895, -0.0627434, 0.18823},
             {0., 0., -0.367052, -0.100618, 1.10116, -0.173587, -0.0475844, -0.0308044, -0.0714449, 0.214335},
             {0., 0., -0.0000513431, -0.0000140744, 0.000154029, -0.0000947733, -0.0000259796, -0.0000117137,
              -9.99369e-6, 0.0000299811}}};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 3; idx2 < 11; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), xi * ffEval[idx1][idx2 - 1]);
                EXPECT_NEAR(comRe, 1., 2e-3);
            }
        }
    }

    TEST(FFBtoDBLPRXPVarTest, evalLOShift) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoDBLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoDBLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        set.changeSetting<bool>("BtoDBLPRXPVar", "WithAsOverM", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        double w = (2 * 5280 * 2381.3371929065643) / (2. * 5280. * 1869);

        // Compare to direct evaluation
        auto vLinXi = ff.XiLinearCoeffs(w);

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            double comRe = compareVals(t.element({idx1, 1}).real(), vLinXi[0] * t.element({idx1, 0}).real());
            EXPECT_NEAR(comRe, 1., 2e-3);
            comRe = compareVals(t.element({idx1, 2}).real(), vLinXi[1] * t.element({idx1, 0}).real());
            EXPECT_NEAR(comRe, 1., 2e-3);
        }
    }

    TEST(FFBtoDBLPRXPVarTest, evalAsOM) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXPVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXPVar", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXPVar", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXPVar", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXPVar", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXPVar", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXPVar", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXPVar", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXPVar", "la2", 0.13);
        set.changeSetting<double>("BtoDBLPRXPVar", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoDBLPRXPVar", "cSt", 1.55);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});
        auto t = ff.getTensor();

        set.changeSetting<bool>("BtoDBLPRXPVar", "WithAsOverM", false);

        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});
        auto t2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2381.3371929065643) / (2. * 5280. * 1869));
        // Compare to direct evaluation
        array<array<double, 10>, 4> ffEval = {
            {{0., 0., -421.757, -115.614, 1265.27, 25.5863, 7.01382, 0., 0., 0.},
             {0., 0., -0.0852484, -0.0233686, 0.255745, -0.00333704, -0.000914761, 0., 0., 0.},
             {0., 0., -0.0989362, -0.0271208, 0.296809, -0.0114073, -0.00312702, 0., 0., 0.},
             {0., 0., -0.0000197618, -5.41718e-6, 0.0000592854, -6.57807e-6, -1.8032e-6, 0., 0., 0.}}};


        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 3; idx2 < 11; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real() - t2.element({idx1, idx2}).real(),
                                           xi * ffEval[idx1][idx2 - 1]);
                EXPECT_NEAR(comRe, 1., 2e-3);
            }
        }
    }

    TEST(FFBtoDBLPRXPVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBLPRXPVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DPLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_TRUE(std::isfinite(t.element({idx1, idx2}).real()))
                    << "element (" << idx1 << "," << idx2 << ") is not finite";
            }
        }
    }

    TEST(FFBtoDBLPRXPVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBLPRXPVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DPLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1.;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDBLPRXPVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBLPRXPVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

} // namespace Hammer
