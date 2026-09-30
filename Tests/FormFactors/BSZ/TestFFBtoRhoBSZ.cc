///
/// @file  TestFFBtoRhoBSZ.cc
/// @brief Tests for FFBtoRhoBSZ
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                   \
    FRIEND_TEST(FFBtoRhoBSZTest, eval);          \
    FRIEND_TEST(FFBtoRhoBSZTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoRhoBSZTest, uninitialized); \
    FRIEND_TEST(FFBtoRhoBSZTest, clone);         \
    FRIEND_TEST(FFBtoRhoBSZTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BSZ/FFBtoRhoBSZ.hh"
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

    TEST(FFBtoRhoBSZTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{17838.5064935065, 0, 0, 17039.1875956242};
        FourMomentum pRhomes{770., 0, 0, 0};
        FourMomentum kNuMu{198.6188302152, 31.769954441031, 97.777865759968, 169.9399855193};
        FourMomentum pMu{16869.8876632913, -31.769954441031, -97.777865759968, 16869.2476101049};

        // Evaluate the FF class
        FFBtoRhoBSZ ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 521), {Particle(pRhomes, 113), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Ap,V,A0,A1,A12,T1,T2,T23
        Tensor ffEval{"ffEval",
                      MD::makeVector({8}, {FF_BRHO},
                                     {-0.12031819247473, 0.34558307298139, 0.37540838626044, 0.26721909723634,
                                      0.30077407846314, 0.28712896725280, 0.27709651260590, 0.7574294922472})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; idx1++) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoRhoBSZTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoRhoBSZ ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BPLUS);
        const double Mc = pdg.getMass(PID::RHO0);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoRhoBSZTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoRhoBSZ ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BPLUS);
        const double Mc = pdg.getMass(PID::RHO0);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoRhoBSZTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoRhoBSZ ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{17838.5064935065, 0, 0, 17039.1875956242};
        FourMomentum pRhomes{770., 0, 0, 0};
        FourMomentum kNuMu{198.6188302152, 31.769954441031, 97.777865759968, 169.9399855193};
        FourMomentum pMu{16869.8876632913, -31.769954441031, -97.777865759968, 16869.2476101049};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 521), {Particle(pRhomes, 113), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBmes, 521), {Particle(pRhomes, 113), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFBtoRhoBSZTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoRhoBSZ ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Straub:2015ica"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Straub:2015ica"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Straub:2015ica"));

        const string tmpfile = testing::TempDir() + "hammer_bsz_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Straub:2015ica"), string::npos);
    }

} // namespace Hammer
