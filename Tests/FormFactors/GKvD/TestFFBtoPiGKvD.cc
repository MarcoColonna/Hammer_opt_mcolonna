///
/// @file  TestFFBtoPiGKvD.cc
/// @brief Tests for FFBtoPiGKvD
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                   \
    FRIEND_TEST(FFBtoPiGKvDTest, eval);          \
    FRIEND_TEST(FFBtoPiGKvDTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoPiGKvDTest, uninitialized); \
    FRIEND_TEST(FFBtoPiGKvDTest, clone);         \
    FRIEND_TEST(FFBtoPiGKvDTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/GKvD/FFBtoPiGKvD.hh"
#include "Hammer/Math/Constants.hh"
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

    TEST(FFBtoPiGKvDTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs (q^2 = 19.8958 GeV^2)
        FourMomentum pBmes{28783.8884892086, 0, 0, 28295.4737822005};
        FourMomentum pPimes{139., 0, 0, 0};
        FourMomentum kNuMu{483.601218826, 143.20931144041, 440.75294019035, 138.1957947273};
        FourMomentum pMu{28161.2872703827, -143.20931144041, -440.75294019035, 28157.2779874733};

        // Evaluate the FF class
        FFBtoPiGKvD ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pPimes, -211), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Fs, Fz, Fp, Ft
        Tensor ffEval{
            "ffEval",
            MD::makeVector({4}, {FF_BPI}, {3109.8859405624, 0.5263315492142, 1.3230756246809, 0.00016462135373268})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoPiGKvDTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoPiGKvD ff;
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

    TEST(FFBtoPiGKvDTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoPiGKvD ff;
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

    TEST(FFBtoPiGKvDTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiGKvD ff;
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

    TEST(FFBtoPiGKvDTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoPiGKvD ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Gubernari:2018wyi"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Gubernari:2018wyi"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Gubernari:2018wyi"));

        const string tmpfile = testing::TempDir() + "hammer_gkvd_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Gubernari:2018wyi"), string::npos);
    }


} // namespace Hammer
