///
/// @file  TestFFBtoOmegaBSZ.cc
/// @brief Tests for FFBtoOmegaBSZ
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                     \
    FRIEND_TEST(FFBtoOmegaBSZTest, eval);          \
    FRIEND_TEST(FFBtoOmegaBSZTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoOmegaBSZTest, uninitialized); \
    FRIEND_TEST(FFBtoOmegaBSZTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BSZ/FFBtoOmegaBSZ.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBtoOmegaBSZTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{17576.6777493606, 0, 0, 16764.8799788388};
        FourMomentum pOmegames{782., 0, 0, 0};
        FourMomentum kNuMu{195.8915304338, 31.769954441031, 97.777865759968, 166.7442672416};
        FourMomentum pMu{16598.7862189268, -31.769954441031, -97.777865759968, 16598.1357115972};

        // Evaluate the FF class
        FFBtoOmegaBSZ ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 521), {Particle(pOmegames, 223), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        // Compare to direct evaluation
        // Ap,V,A0,A1,A12,T1,T2,T23
        Tensor ffEval{"ffEval",
                      MD::makeVector({8}, {FF_BOMEGA},
                                     {-0.11262003760786, 0.32133945171369, 0.34599698254843, 0.24779174538358,
                                      0.27431431527030, 0.26569423684374, 0.25639515670707, 0.69258666245592})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; idx1++) {
            EXPECT_NEAR(compareVals(t.element({idx1}).real(), ffEval.element({idx1}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(t.element({idx1}).imag(), ffEval.element({idx1}).imag());
        };
    }

    TEST(FFBtoOmegaBSZTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoOmegaBSZ ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BPLUS);
        const double Mc = pdg.getMass(PID::OMEGA);
        const double Sqq = Mb * Mb + Mc * Mc - (2. * Mb * Mc);

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});

        auto& t = ff.getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_TRUE(std::isfinite(t.element({idx}).real()));
        }
    }

    TEST(FFBtoOmegaBSZTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoOmegaBSZ ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BPLUS);
        const double Mc = pdg.getMass(PID::OMEGA);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoOmegaBSZTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoOmegaBSZ ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{17576.6777493606, 0, 0, 16764.8799788388};
        FourMomentum pOmegames{782., 0, 0, 0};
        FourMomentum kNuMu{195.8915304338, 31.769954441031, 97.777865759968, 166.7442672416};
        FourMomentum pMu{16598.7862189268, -31.769954441031, -97.777865759968, 16598.1357115972};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 521), {Particle(pOmegames, 223), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        ff.eval(Particle(pBmes, 521), {Particle(pOmegames, 223), Particle(kNuMu, 14), Particle(pMu, -13)}, {});

        auto& tOrig = ff.getTensor();
        auto& tCloned = cloned->getTensor();
        for (IndexType idx = 0; idx < 8; ++idx) {
            EXPECT_NEAR(compareVals(tOrig.element({idx}).real(), tCloned.element({idx}).real()), 1., 1.e-4);
            EXPECT_DOUBLE_EQ(tOrig.element({idx}).imag(), tCloned.element({idx}).imag());
        }
    }

} // namespace Hammer
