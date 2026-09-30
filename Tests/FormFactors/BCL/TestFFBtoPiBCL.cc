///
/// @file  TestFFBtoPiBCL.cc
/// @brief Tests for FFBtoPiBCL
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                  \
    FRIEND_TEST(FFBtoPiBCLTest, eval);          \
    FRIEND_TEST(FFBtoPiBCLTest, evalq2cons);    \
    FRIEND_TEST(FFBtoPiBCLTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoPiBCLTest, uninitialized); \
    FRIEND_TEST(FFBtoPiBCLTest, clone);         \
    FRIEND_TEST(FFBtoPiBCLTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BCL/FFBtoPiBCL.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoPiBCLTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 = 19.8958 GeV^2)
        FourMomentum pBmes{28783.8884892086, 0, 0, 28295.4737822005};
        FourMomentum pPimes{139., 0, 0, 0};
        FourMomentum kNuMu{483.601218826, 143.20931144041, 440.75294019035, 138.1957947273};
        FourMomentum pMu{28161.2872703827, -143.20931144041, -440.75294019035, 28157.2779874733};

        // Evaluate the FF class
        FFBtoPiBCL ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BPI}, {0., 0.50165560937469, 1.3961877732509, 0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 3; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoPiBCLTest, evalq2cons) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 = 19.8958 GeV^2)
        FourMomentum pBmes{28783.8884892086, 0, 0, 28295.4737822005};
        FourMomentum pPimes{139., 0, 0, 0};
        FourMomentum kNuMu{483.601218826, 143.20931144041, 440.75294019035, 138.1957947273};
        FourMomentum pMu{28161.2872703827, -143.20931144041, -440.75294019035, 28157.2779874733};

        // Evaluate the FF class
        FFBtoPiBCL ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        vector<double> apvec = {0.415, -0.488, -0.31};
        vector<double> a0vec = {0.500, -1.424};
        set.changeSetting<vector<double>>("BtoPiBCL", "ap", apvec);
        set.changeSetting<vector<double>>("BtoPiBCL", "a0", a0vec);
        set.changeSetting<bool>("BtoPiBCL", "q2cons", true);
        ff.calcUnits();

        ff.eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BPI}, {0., 0.49302579294957, 1.3829059390230, 0.})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 3; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoPiBCLTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoPiBCL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::PIMINUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoPiBCLTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoPiBCL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::PIMINUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoPiBCLTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiBCL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{28783.8884892086, 0, 0, 28295.4737822005};
        FourMomentum pPimes{139., 0, 0, 0};
        FourMomentum kNuMu{483.601218826, 143.20931144041, 440.75294019035, 138.1957947273};
        FourMomentum pMu{28161.2872703827, -143.20931144041, -440.75294019035, 28157.2779874733};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFBtoPiBCLTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiBCL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Aoki:2019cca"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Aoki:2019cca"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Aoki:2019cca"));

        const string tmpfile = testing::TempDir() + "hammer_bcl_btopi_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Aoki:2019cca"), string::npos);
    }

} // namespace Hammer
