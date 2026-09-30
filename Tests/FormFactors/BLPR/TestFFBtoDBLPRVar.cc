///
/// @file  TestFFBtoDBLPRVar.cc
/// @brief Tests for FFBtoDBLPRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFBtoDBLPRVarTest, eval);          \
    FRIEND_TEST(FFBtoDBLPRVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoDBLPRVarTest, uninitialized); \
    FRIEND_TEST(FFBtoDBLPRVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLPR/FFBtoDBLPRVar.hh"
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

    TEST(FFBtoDBLPRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2945.70312500000, 0, 0, -2280.11883476164};
        FourMomentum pTau{2284.58811187298, 24.8428481565870, 43.029075211928, 2281.63300565477};
        FourMomentum kNuTau{49.70876312702, -24.8428481565870, -43.029075211928, -1.51417089313};

        // Evaluate the FF class
        FFBtoDBLPRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        Tensor ffEval{
            "ffEval",
            MD::makeVector(
                {4, 8}, {FF_BD, FF_BD_VAR},
                {4466.671011134957,      -2470.24429911416,      -2684.96443646855,      -1295.05495659087,
                 8054.89330940564,       -503.739261988966,      -303.525086179545,      -2.79690796197956,
                 0.6703585240318198,     -0.370734562322605,     -0.376654896744685,     -0.179119228472118,
                 1.12996469023405,       -0.0833100738272369,    -0.0500200793272945,    -0.000419760284241094,
                 0.6739387878244352,     -0.372714588655636,     -0.378513905576002,     -0.179987424080559,
                 1.13554171672801,       -0.0851350920385681,    -0.0510869328722607,    -0.000422002148099567,
                 0.00010824070211699414, -0.0000598613546128796, -0.0000532446680547304, -0.0000245337491907538,
                 0.000159734004164191,   -0.0000441403198753062, -0.0000258595138797868, -6.77773851726666e-8})};


        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 2e-3);
            }
        }
    }

    TEST(FFBtoDBLPRVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBLPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DMINUS);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                EXPECT_TRUE(std::isfinite(t.element({idx1, idx2}).real()))
                    << "element (" << idx1 << "," << idx2 << ") is not finite";
            }
        }
    }

    TEST(FFBtoDBLPRVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoDBLPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DMINUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoDBLPRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoDBLPRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2945.70312500000, 0, 0, -2280.11883476164};
        FourMomentum pTau{2284.58811187298, 24.8428481565870, 43.029075211928, 2281.63300565477};
        FourMomentum kNuTau{49.70876312702, -24.8428481565870, -43.029075211928, -1.51417089313};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pDmes, -411), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx1 = 0; idx1 < 4; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }


} // namespace Hammer
