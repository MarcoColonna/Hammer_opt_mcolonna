///
/// @file  TestFFBctoJpsiKiselev.cc
/// @brief Tests for FFBctoJpsiKiselev
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                                      \
    FRIEND_TEST(FFBctoJpsiKiselevTest, eval);                       \
    FRIEND_TEST(FFBctoJpsiKiselevTest, zeroRecoil);                 \
    FRIEND_TEST(FFBctoJpsiKiselevTest, uninitialized);              \
    FRIEND_TEST(FFBctoJpsiKiselevTest, addRefs);                    \
    FRIEND_TEST(FFBctoJpsiKiselevTest, evalAtPSPointDirect);        \
    FRIEND_TEST(FFBctoJpsiKiselevTest, evalAtPSPointUninitialized); \
    FRIEND_TEST(FFBctoJpsiKiselevTest, evalAtPSPointSecondPoint);   \
    FRIEND_TEST(FFBctoJpsiKiselevTest, clone);                      \
    FRIEND_TEST(FFBctoJpsiKiselevTest, getFFPSIntegrandNonEmpty)


#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/Kiselev/FFBctoJpsiKiselev.hh"
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

    TEST(FFBctoJpsiKiselevTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBcmes{6275., 0., 0., 0.};
        FourMomentum pJpsimes{3821.5809561753, -857.35310207508, 0., -2069.83348677229};
        FourMomentum kNuTau{131.0024243058, -0.827618332, 98.867858366606, -85.9424039919};
        FourMomentum pTau{2322.4166195189, 858.1807204071, -98.867858366606, 2155.7758907642};

        // Evaluate the FF class
        FFBctoJpsiKiselev ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443), Particle(kNuTau, 12), Particle(pTau, -11)}, {});
        // ff.evalAtPSPoint({Sqq}, {pBmes.mass(), pJpsimes.mass()});

        // Compare to direct evaluation
        // Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt
        Tensor ffEval{"ffEval", MD::makeVector({8}, {FF_BCJPSI},
                                               {0., 6206.49, 0.000115714, 0.000126234, -0.0000778442, 0., 0., 0.})};

        auto& t = ff.getTensor();
        EXPECT_TRUE(t.isSameLabelShape(ffEval));
        for (IndexType idx1 = 1; idx1 < 5; idx1++) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        };
    }

    TEST(FFBctoJpsiKiselevTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBctoJpsiKiselev ff;
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

    TEST(FFBctoJpsiKiselevTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBctoJpsiKiselev ff;
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

    TEST(FFBctoJpsiKiselevTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        FFBctoJpsiKiselev ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Kiselev:2002vz"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Kiselev:2002vz"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Kiselev:2002vz"));
    }

    TEST(FFBctoJpsiKiselevTest, evalAtPSPointDirect) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiKiselev ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 6.2749;
        const double Mc = 3.09688;
        const double t = 3.0; // GeV^2

        static_cast<FormFactorBase*>(&ff)->evalAtPSPoint({t}, {Mb, Mc});

        const double sumMasses = Mb + Mc;
        const double diffMasses = Mb - Mc;
        const double Mpole2 = 4.5 * 4.5;
        const double Den = 1.0 / (1.0 - t / Mpole2);
        const double FV = 0.11 * Den;
        const double FAp = -0.074 * Den;
        const double FA0 = 5.9 * Den;
        const double FAm = 0.12 * Den;

        const double Vf = sumMasses * FV;
        const double A2f = -sumMasses * FAp;
        const double A1f = FA0 / sumMasses;
        const double A0f = (t * FAm + sumMasses * A1f - diffMasses * A2f) / (sumMasses - diffMasses);
        const double rC = Mc / Mb;

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

    TEST(FFBctoJpsiKiselevTest, evalAtPSPointSecondPoint) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiKiselev ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 6.2749;
        const double Mc = 3.09688;
        const double t = 8.0; // GeV^2

        static_cast<FormFactorBase*>(&ff)->evalAtPSPoint({t}, {Mb, Mc});

        const double sumMasses = Mb + Mc;
        const double diffMasses = Mb - Mc;
        const double Mpole2 = 4.5 * 4.5;
        const double Den = 1.0 / (1.0 - t / Mpole2);
        const double FV = 0.11 * Den;
        const double FAp = -0.074 * Den;
        const double FA0 = 5.9 * Den;
        const double FAm = 0.12 * Den;

        const double Vf = sumMasses * FV;
        const double A2f = -sumMasses * FAp;
        const double A1f = FA0 / sumMasses;
        const double A0f = (t * FAm + sumMasses * A1f - diffMasses * A2f) / (sumMasses - diffMasses);
        const double rC = Mc / Mb;

        const double expFf = A1f * Mb * (1.0 + rC);
        const double expFg = Vf / (Mb * (1.0 + rC));
        const double expFm = (A2f * (1.0 - rC) + 2.0 * A0f * rC - A1f * (1.0 + rC)) * Mb / t;
        const double expFp = -A2f / (Mb * (1.0 + rC));

        auto& ten = ff.getTensor();
        EXPECT_NEAR(compareVals(ten.element({1}).real(), expFf), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({2}).real(), expFg), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({3}).real(), expFm), 1., 1e-10);
        EXPECT_NEAR(compareVals(ten.element({4}).real(), expFp), 1., 1e-10);
    }

    TEST(FFBctoJpsiKiselevTest, evalAtPSPointUninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiKiselev ff;
        ff.setSettingsHandler(set);
        // Deliberately skip initSettings() — _initialized remains false.
        EXPECT_NO_THROW(static_cast<FormFactorBase*>(&ff)->evalAtPSPoint({3.0}, {6.2749, 3.09688}));
    }

    TEST(FFBctoJpsiKiselevTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        FFBctoJpsiKiselev ff;
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

    TEST(FFBctoJpsiKiselevTest, getFFPSIntegrandNonEmpty) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");
        FFBctoJpsiKiselev ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        const double Mb = 6.2749;
        const double Mc = 3.09688;

        EvaluationGrid grid{{1.0, Mb, Mc}, {3.0, Mb, Mc}, {6.0, Mb, Mc}};

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
