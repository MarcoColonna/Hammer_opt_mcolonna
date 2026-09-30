///
/// @file  TestFFBctoJpsiEFG.cc
/// @brief Tests for FFBctoJpsiEFG
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                                  \
    FRIEND_TEST(FFBctoJpsiEFGTest, eval);                       \
    FRIEND_TEST(FFBctoJpsiEFGTest, zeroRecoil);                 \
    FRIEND_TEST(FFBctoJpsiEFGTest, uninitialized);              \
    FRIEND_TEST(FFBctoJpsiEFGTest, addRefs);                    \
    FRIEND_TEST(FFBctoJpsiEFGTest, evalAtPSPointDirect);        \
    FRIEND_TEST(FFBctoJpsiEFGTest, evalAtPSPointUninitialized); \
    FRIEND_TEST(FFBctoJpsiEFGTest, evalAtPSPointTZero);         \
    FRIEND_TEST(FFBctoJpsiEFGTest, clone);                      \
    FRIEND_TEST(FFBctoJpsiEFGTest, getFFPSIntegrandNonEmpty)


#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/EFG/FFBctoJpsiEFG.hh"
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

    TEST(FFBctoJpsiEFGTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBcmes{6275., 0., 0., 0.};
        FourMomentum pJpsimes{3821.5809561753, -857.35310207508, 0., -2069.83348677229};
        FourMomentum kNuTau{131.0024243058, -0.827618332, 98.867858366606, -85.9424039919};
        FourMomentum pTau{2322.4166195189, 858.1807204071, -98.867858366606, 2155.7758907642};

        // Evaluate the FF class
        FFBctoJpsiEFG ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443), Particle(kNuTau, 12), Particle(pTau, -11)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pJpsimes.mass()});

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BCJPSI},
                                               {0., 4838.60027371874, 0.0000557416471047622, 0.000337714145792924,
                                                -0.0000821548486986795, 0., 0., 0.})};

        auto& t = ff.getTensor();
        EXPECT_TRUE(t.isSameLabelShape(ffEval));
        for (IndexType idx1 = 1; idx1 < 5; idx1++) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        };
    }

    TEST(FFBctoJpsiEFGTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBctoJpsiEFG ff;
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

    TEST(FFBctoJpsiEFGTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBctoJpsiEFG ff;
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

    TEST(FFBctoJpsiEFGTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        FFBctoJpsiEFG ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Ebert:2003cn"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Ebert:2003cn"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Ebert:2003cn"));
    }

    TEST(FFBctoJpsiEFGTest, evalAtPSPointDirect) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiEFG ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double t = 5.0; // GeV^2
        const double Mb = 6.2749;
        const double Mc = 3.09688;

        static_cast<FormFactorBase*>(&ff)->evalAtPSPoint({t}, {Mb, Mc});

        const double rC = Mc / Mb;
        const double Vf = (0.49077824756158533 - 0.0012925655191347828 * t) / (1.0 - 0.06292520325875656 * t);
        const double A0f = (0.4160345034630221 - 0.0024720095310225023 * t) / (1.0 - 0.061603451915567785 * t);
        const double A1f = (0.4970212860605933 - 0.0067519730024654745 * t) / (1.0 - 0.050487026667172176 * t);
        const double A2f = (0.7315284919705497 + 0.0014263826220727142 * t - 0.0006946090066269195 * t * t) /
                           (1.0 - 0.04885587273651653 * t);

        const double expFf = A1f * Mb * (1.0 + rC);
        const double expFg = Vf / (Mb * (1.0 + rC));
        const double expFm = (A2f * (1.0 - rC) + 2.0 * A0f * rC - A1f * (1.0 + rC)) * Mb / t;
        const double expFp = -A2f / (Mb * (1.0 + rC));

        auto& ten = ff.getTensor();
        EXPECT_NEAR(compareVals(ten.element({1}).real(), expFf), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({2}).real(), expFg), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({3}).real(), expFm), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({4}).real(), expFp), 1., 1e-10);

        EXPECT_DOUBLE_EQ(ten.element({0}).real(), 0.);
        EXPECT_DOUBLE_EQ(ten.element({5}).real(), 0.);
        EXPECT_DOUBLE_EQ(ten.element({6}).real(), 0.);
        EXPECT_DOUBLE_EQ(ten.element({7}).real(), 0.);
    }

    TEST(FFBctoJpsiEFGTest, evalAtPSPointTZero) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiEFG ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 6.2749;
        const double Mc = 3.09688;
        const double Sqq = 0.0;

        static_cast<FormFactorBase*>(&ff)->evalAtPSPoint({Sqq}, {Mb, Mc});

        const double rC = Mc / Mb;
        const double Vf = 0.49077824756158533;
        const double A1f = 0.4970212860605933;
        const double A2f = 0.7315284919705497;

        const double expFf = A1f * Mb * (1.0 + rC);
        const double expFg = Vf / (Mb * (1.0 + rC));
        const double expFp = -A2f / (Mb * (1.0 + rC));

        auto& ten = ff.getTensor();
        EXPECT_NEAR(compareVals(ten.element({1}).real(), expFf), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({2}).real(), expFg), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({4}).real(), expFp), 1., 1e-10);
    }

    TEST(FFBctoJpsiEFGTest, evalAtPSPointUninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiEFG ff;
        ff.setSettingsHandler(set);

        EXPECT_NO_THROW(static_cast<FormFactorBase*>(&ff)->evalAtPSPoint({5.0}, {6.2749, 3.09688}));
    }

    TEST(FFBctoJpsiEFGTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        FFBctoJpsiEFG ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("myLabel");
        ASSERT_NE(cloned, nullptr);

        EXPECT_NE(cloned->group().find("myLabel"), string::npos);

        FourMomentum pBcmes{6275., 0., 0., 0.};
        FourMomentum pJpsimes{3821.5809561753, -857.35310207508, 0., -2069.83348677229};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443)}, {});

        ff.eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 1; idx < 5; ++idx) {
            EXPECT_NEAR(tOrig.element({idx}).real(), tCloned.element({idx}).real(), 1e-12);
        }
    }

    TEST(FFBctoJpsiEFGTest, getFFPSIntegrandNonEmpty) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiEFG ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 6.2749;
        const double Mc = 3.09688;

        EvaluationGrid grid{{0.5, Mb, Mc}, {3.0, Mb, Mc}, {7.0, Mb, Mc}};

        Tensor integrand = ff.getFFPSIntegrand(grid);

        EXPECT_GT(integrand.rank(), 0u);

        bool anyNonZero = false;
        for (IndexType i = 0; i < static_cast<IndexType>(grid.size()); ++i) {
            for (IndexType j = 1; j < 5; ++j) {
                if (!isZero(integrand.element({j, j, i}).real())) {
                    anyNonZero = true;
                }
            }
        }
        EXPECT_TRUE(anyNonZero);
    }

} // namespace Hammer
