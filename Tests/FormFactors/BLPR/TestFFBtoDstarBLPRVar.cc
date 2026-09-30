///
/// @file  TestFFBtoDstarBLPRVar.cc
/// @brief Tests for FFBtoDstarBLPRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                         \
    FRIEND_TEST(FFBtoDstarBLPRVarTest, eval);          \
    FRIEND_TEST(FFBtoDstarBLPRVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDstarBLPRVarTest, uninitialized); \
    FRIEND_TEST(FFBtoDstarBLPRVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLPR/FFBtoDstarBLPRVar.hh"
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

    TEST(FFBtoDstarBLPRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2998.91098484848, 0., 0., -2225.61611583061};
        FourMomentum pTau{2231.33061982896, 24.842848156587, 43.029075211928, 2228.30488635828};
        FourMomentum kNuTau{49.75839532255, -24.842848156587, -43.029075211928, -2.68877052767};

        // Evaluate the FF class
        FFBtoDstarBLPRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        Tensor ffEval{
            "ffEval",
            MD::makeVector(
                {8, 8}, {FF_BD, FF_BD_VAR},
                {-0.45031812876186633,   0.213641396192117,      -0.104749873769849,    -0.0482927137349265,
                 0.0270901546545336,     -0.105331922997079,     -0.0519717407132333,   -0.0000358244106630322,
                 4534.436859561684,      -2151.2423323602,       -514.323205500589,     -285.707395417011,
                 -346.214655288905,      -105.425934127357,      -50.3696811550195,     0.360730598675119,
                 0.00010688127403113114, -0.0000507069628166762, -9.67176811555882e-6,  -5.52835913713112e-6,
                 -6.70104722733027e-6,   -0.0000100954652565242, -4.93158111673046e-6,  8.50278593848067e-9,
                 0.00010031462234340629, -0.0000475915904937139, -0.000109918290211195, -0.0000548018971496432,
                 -6.65274178487814e-6,   -0.0000191162126968543, -9.37191979752968e-6,  7.98038541379168e-9,
                 -0.0000792681778066827, 0.0000376066675947219,  0.0000546988013669274, 0.0000274825534912616,
                 6.49792044460652e-6,    0.0000302059389079022,  0.000014834974797635,  -6.30606580744283e-9,
                 -3.2486347178931773e-9, 1.54122788428956e-9,    -1.37569946876016e-8,  -6.74497895406402e-9,
                 2.38975272141354e-11,   -2.7586383609572e-9,    -1.35831192912498e-9,  -2.58440459743369e-13,
                 0.4069959622372874,     -0.193088352574207,     -0.0314857026049396,   -0.018422518381524,
                 -0.0223353618321834,    0.0735959417200616,     0.0363434548927844,    0.0000323779780518116,
                 -0.6834765758059944,    0.324257187516021,      0.0707418219210209,    0.0397279121670139,
                 0.0481467372114464,     -0.0330121713888342,    -0.0164678459627562,   -0.000054372995370091})};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 2e-3);
            }
        }
    }

    TEST(FFBtoDstarBLPRVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBLPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSTARMINUS);
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

    TEST(FFBtoDstarBLPRVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDstarBLPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSTARMINUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDstarBLPRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDstarBLPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2998.91098484848, 0., 0., -2225.61611583061};
        FourMomentum pTau{2231.33061982896, 24.842848156587, 43.029075211928, 2228.30488635828};
        FourMomentum kNuTau{49.75839532255, -24.842848156587, -43.029075211928, -2.68877052767};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

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
