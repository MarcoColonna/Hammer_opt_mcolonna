///
/// @file  TestFFBtoDISGW2.cc
/// @brief Tests for FFBtoDISGW2
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                   \
    FRIEND_TEST(FFBtoDISGW2Test, eval);          \
    FRIEND_TEST(FFBtoDISGW2Test, uninitialized); \
    FRIEND_TEST(FFBtoDISGW2Test, bsMode);        \
    FRIEND_TEST(FFBtoDISGW2Test, clone);         \
    FRIEND_TEST(FFBtoDISGW2Test, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/ISGW2/FFBtoDISGW2.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoDISGW2Test, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2381.3371929065643, 0., 0., -1475.671313782346};
        FourMomentum pTau{2034.777874233763, 553.6990600712329, -92.91177555543321, 819.1030156888944};
        FourMomentum kNuTau{863.8849328596731, -553.6990600712329, 92.91177555543321, 656.5682980934517};


        // Evaluate the FF class
        FFBtoDISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        //        TestableTensor FFres;
        //        dynamic_cast<Tensor&>(FFres)=ff.evalAtPSPoint({7.6762},{5.280,1.86961});
        // Switched to eval on momenta (q^2 = 6.22464 GeV^2)
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 16), Particle(pTau, -15)}, {});

        // Compare to direct evaluation
        // The calculated form factors at t=7.6762 are f+=1.03227 f0=0.862739
        // Fs, Fz, Fp, Ft
        Tensor ffEval{"ffEval", MD::makeVector({4}, {FF_BD}, {0, 0.856145172884, 0.98763713683, 0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 1; idx1 < 3; ++idx1) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        }
    }

    TEST(FFBtoDISGW2Test, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5280.; // B0 (MeV)
        const double Mc = 1869.; // D- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDISGW2Test, bsMode) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDISGW2 ff;
        ff.setSettingsHandler(set);
        ff.setSignatureIndex(2); // Bs -> Ds- signature
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 5367.; // Bs (MeV)
        const double Mc = 1968.; // Ds- (MeV)
        const double Sqq = 0.5 * (Mb - Mc) * (Mb - Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 4; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoDISGW2Test, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDISGW2 ff;
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

    TEST(FFBtoDISGW2Test, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDISGW2 ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Isgur:1988gb"));
        EXPECT_FALSE(set.checkReference("Scora:1995ty"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Isgur:1988gb"));
        EXPECT_TRUE(set.checkReference("Scora:1995ty"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Isgur:1988gb"));
        EXPECT_TRUE(set.checkReference("Scora:1995ty"));

        const string tmpfile = testing::TempDir() + "hammer_isgw2_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Isgur:1988gb"), string::npos);
        EXPECT_NE(content.find("Scora:1995ty"), string::npos);
    }

} // namespace Hammer
