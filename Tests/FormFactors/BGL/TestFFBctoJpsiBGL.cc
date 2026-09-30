///
/// @file  TestFFBctoJpsiBGL.cc
/// @brief Tests for FFBctoJpsiBGL
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFBctoJpsiBGLTest, eval);          \
    FRIEND_TEST(FFBctoJpsiBGLTest, zeroRecoil);    \
    FRIEND_TEST(FFBctoJpsiBGLTest, uninitialized); \
    FRIEND_TEST(FFBctoJpsiBGLTest, clone);         \
    FRIEND_TEST(FFBctoJpsiBGLTest, addRefs)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BGL/FFBctoJpsiBGL.hh"
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

    TEST(FFBctoJpsiBGLTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBcmes{6275., 0., 0., 0.};
        FourMomentum pJpsimes{3821.5809561753, -857.35310207508, 0., -2069.83348677229};
        FourMomentum kNuTau{131.0024243058, -0.827618332, 98.867858366606, -85.9424039919};
        FourMomentum pTau{2322.4166195189, 858.1807204071, -98.867858366606, 2155.7758907642};

        // Evaluate the FF class
        FFBctoJpsiBGL ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443), Particle(kNuTau, 12), Particle(pTau, -11)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pJpsimes.mass()});

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BCJPSI},
                                               {0., 4427.03290506073, 0.0000826899544190131, 0.000250864518586240,
                                                -0.0000481165205547890, 0., 0., 0.})};

        auto& t = ff.getTensor();
        EXPECT_TRUE(t.isSameLabelShape(ffEval));
        for (IndexType idx1 = 1; idx1 < 5; idx1++) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        };
    }

    TEST(FFBctoJpsiBGLTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBctoJpsiBGL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BCPLUS);
        const double Mc = pdg.getMass(PID::JPSI);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBctoJpsiBGLTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBctoJpsiBGL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BCPLUS);
        const double Mc = pdg.getMass(PID::JPSI);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBctoJpsiBGLTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBctoJpsiBGL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBcmes{6275., 0., 0., 0.};
        FourMomentum pJpsimes{3821.5809561753, -857.35310207508, 0., -2069.83348677229};
        FourMomentum kNuTau{131.0024243058, -0.827618332, 98.867858366606, -85.9424039919};
        FourMomentum pTau{2322.4166195189, 858.1807204071, -98.867858366606, 2155.7758907642};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443), Particle(kNuTau, 12), Particle(pTau, -11)}, {});

        ff.eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443), Particle(kNuTau, 12), Particle(pTau, -11)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

    TEST(FFBctoJpsiBGLTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBctoJpsiBGL ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Cohen:2019zev"));
        EXPECT_FALSE(set.checkReference("Harrison:2020gvo"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Cohen:2019zev"));
        EXPECT_TRUE(set.checkReference("Harrison:2020gvo"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Cohen:2019zev"));
        EXPECT_TRUE(set.checkReference("Harrison:2020gvo"));

        const string tmpfile = testing::TempDir() + "hammer_bgl_FFBctoJpsiBGL_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Cohen:2019zev"), string::npos);
        EXPECT_NE(content.find("Harrison:2020gvo"), string::npos);
    }


} // namespace Hammer
