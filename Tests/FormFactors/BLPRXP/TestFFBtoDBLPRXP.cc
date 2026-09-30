///
/// @file  TestFFBtoDBLPRXP.cc
/// @brief Tests for FFBtoDBLPRXP
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                    \
    FRIEND_TEST(FFBtoDBLPRXPTest, AsOverM);       \
    FRIEND_TEST(FFBtoDBLPRXPTest, OneOverMbMc);   \
    FRIEND_TEST(FFBtoDBLPRXPTest, OneOverMb2);    \
    FRIEND_TEST(FFBtoDBLPRXPTest, eval);          \
    FRIEND_TEST(FFBtoDBLPRXPTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDBLPRXPTest, uninitialized); \
    FRIEND_TEST(FFBtoDBLPRXPTest, clone);         \
    FRIEND_TEST(FFBtoDBLPRXPTest, addRefs)

#define FFBLPRXPBASE_FRIENDS                      \
    FRIEND_TEST(FFBtoDBLPRXPTest, AsOverM);       \
    FRIEND_TEST(FFBtoDBLPRXPTest, OneOverMbMc);   \
    FRIEND_TEST(FFBtoDBLPRXPTest, OneOverMb2);    \
    FRIEND_TEST(FFBtoDBLPRXPTest, eval);          \
    FRIEND_TEST(FFBtoDBLPRXPTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDBLPRXPTest, uninitialized); \
    FRIEND_TEST(FFBtoDBLPRXPTest, clone);         \
    FRIEND_TEST(FFBtoDBLPRXPTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLPRXP/FFBLPRXPBase.hh"
#include "Hammer/FormFactors/BLPRXP/FFBtoDBLPRXP.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/SettingsConsumer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoDBLPRXPTest, AsOverM) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXP", "la2", 0.13);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});
        auto tDiff = ff.getTensor();

        set.changeSetting<bool>("BtoDBLPRXP", "WithAsOverM", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto tDiff2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2381.3371929065643) / (2. * 5280. * 1869));
        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{
            "ffEval",
            MD::makeVector({4}, {FF_BD},
                           {xi * 7143.91 * 0.0183884, xi * (0.999287 * 0.0250252 + (-0.252456) * (-0.00147448)),
                            xi * (1.13787 * 0.0250252 + (-0.542913) * (-0.00147448)), xi * 0.000159165 * 0.0667049})};

        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            double comRe = compareVals(tDiff.element({idx1}).real() - tDiff2.element({idx1}).real(),
                                       ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDBLPRXPTest, OneOverMbMc) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXP", "la2", 0.13);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});
        auto tDiff = ff.getTensor();

        set.changeSetting<bool>("BtoDBLPRXP", "With1OverMbMc", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto tDiff2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2381.3371929065643) / (2. * 5280. * 1869));
        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{
            "ffEval",
            MD::makeVector({4}, {FF_BD},
                           {xi * 7143.91 * 0.00522267, xi * (0.999287 * (-0.00522267) + (-0.252456) * 0.),
                            xi * (1.13787 * (-0.00522267) + (-0.542913) * 0.), xi * 0.000159165 * 0.00522267})};

        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            double comRe = compareVals(tDiff.element({idx1}).real() - tDiff2.element({idx1}).real(),
                                       ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDBLPRXPTest, OneOverMb2) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXP", "la2", 0.13);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pDstarmes.mass()});
        auto tDiff = ff.getTensor();

        set.changeSetting<bool>("BtoDBLPRXP", "With1OverMb2", false);
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});
        auto tDiff2 = ff.getTensor();

        double xi = ff.Xi((2 * 5280 * 2381.3371929065643) / (2. * 5280. * 1869));
        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BD},
                                               {xi * 7143.91 * 0.000459997,
                                                xi * (0.999287 * 0.000665813 + (-0.252456) * (-0.00170744)),
                                                xi * (1.13787 * 0.000665813 + (-0.542913) * (-0.00170744)),
                                                xi * 0.000159165 * (-0.00104163)})};

        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            double comRe = compareVals(tDiff.element({idx1}).real() - tDiff2.element({idx1}).real(),
                                       ffEval.element({idx1}).real());
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(FFBtoDBLPRXPTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};

        // Evaluate the FF class
        FFBtoDBLPRXP ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        set.changeSetting<double>("BtoDBLPRXP", "chi21", 0.82);
        set.changeSetting<double>("BtoDBLPRXP", "chi2p", -0.69);
        set.changeSetting<double>("BtoDBLPRXP", "chi3p", 0.3);
        set.changeSetting<double>("BtoDBLPRXP", "eta1", 0.34);
        set.changeSetting<double>("BtoDBLPRXP", "phi1p", 0.26);
        set.changeSetting<double>("BtoDBLPRXP", "mb", 4.69);
        set.changeSetting<double>("BtoDBLPRXP", "mc", 4.69 - 3.4);
        set.changeSetting<double>("BtoDBLPRXP", "mDBar", 1.973);
        set.changeSetting<double>("BtoDBLPRXP", "rho1", -0.12);
        set.changeSetting<double>("BtoDBLPRXP", "la2", 0.13);
        set.changeSetting<double>("BtoDBLPRXP", "RhoStSq", 1.47);
        set.changeSetting<double>("BtoDBLPRXP", "cSt", 1.55);
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BD}, {5228.31, 0.783358, 0.8981, 0.0001471})};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoDBLPRXPTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBLPRXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DPLUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoDBLPRXPTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBLPRXP ff;
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

    TEST(FFBtoDBLPRXPTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBLPRXP ff;
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
        cloned->eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFBtoDBLPRXPTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBLPRXP ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Bernlochner:2017jka"));
        EXPECT_FALSE(set.checkReference("Bernlochner:2022ywh"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2017jka"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2022ywh"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Bernlochner:2017jka"));
        EXPECT_TRUE(set.checkReference("Bernlochner:2022ywh"));

        const string tmpfile = testing::TempDir() + "hammer_blprxp_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Bernlochner:2017jka"), string::npos);
        EXPECT_NE(content.find("Bernlochner:2022ywh"), string::npos);
    }

} // namespace Hammer
