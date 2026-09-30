///
/// @file  TestFFBtoD2starBLRVar.cc
/// @brief Tests for FFBtoD2starBLRVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define FORMFACTORBASE_FRIENDS                         \
    FRIEND_TEST(FFBtoD2starBLRVarTest, eval);          \
    FRIEND_TEST(FFBtoD2starBLRVarTest, zeroRecoil);    \
    FRIEND_TEST(FFBtoD2starBLRVarTest, uninitialized); \
    FRIEND_TEST(FFBtoD2starBLRVarTest, clone)

#include "Hammer/FormFactorBase.hh"
#include "Hammer/FormFactors/BLR/FFBtoD2starBLRVar.hh"
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

    TEST(FFBtoD2starBLRVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD2starmes{3189.85994318182, 0, 0, -2029.45447278719};
        FourMomentum kNuTau{50.20101629360, -24.8428481565870, -43.029075211928, -7.17451174540};
        FourMomentum pTau{2039.93904052459, 24.8428481565870, 43.029075211928, 2036.62898453259};

        // Evaluate the FF class
        FFBtoD2starBLRVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBmes, 511), {Particle(pD2starmes, -415), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({8, 8}, {FF_BDSSD2STAR, FF_BDSSD2STAR_VAR},
                                               {0.546490387225108,
                                                0.780700553178726,
                                                0.307620900305189,
                                                0.364527194110234,
                                                0.101473798088436,
                                                -0.281142111174432,
                                                -0.0832642109616469,
                                                0.140571055587216,
                                                -0.688578008829105,
                                                -0.983682869755864,
                                                -0.387602402453807,
                                                -0.0000357979440126166,
                                                0.0532113110391841,
                                                0.645548433310511,
                                                0,
                                                -0.322774216655256,
                                                0.106563686860457,
                                                0.152233838372081,
                                                0.0599849842891214,
                                                -0.281142111174432,
                                                0,
                                                0,
                                                -0.281142111174432,
                                                0,
                                                0.0388647072104949,
                                                0.0555210103007071,
                                                0.0218770476144996,
                                                -0.281021239200277,
                                                -0.179668313085996,
                                                -0.281142111174432,
                                                0.281142111174432,
                                                0.140571055587216,
                                                0.0760885649215852,
                                                0.108697949887979,
                                                0.0428304566580900,
                                                -0.000120871974155222,
                                                0.179668313085996,
                                                0.281142111174432,
                                                0,
                                                -0.140571055587216,
                                                0.379386626757174,
                                                0.541980895367392,
                                                0.213557746695949,
                                                0,
                                                0,
                                                -0.281142111174432,
                                                0,
                                                0.140571055587216,
                                                0.464770127945815,
                                                0.663957325636879,
                                                0.261620347833779,
                                                -0.281021239200277,
                                                0.101473798088436,
                                                0,
                                                0,
                                                0,
                                                -0.114672154190115,
                                                -0.163817363128736,
                                                -0.0645492622313376,
                                                0.281142111174432,
                                                0,
                                                0,
                                                -0.281142111174432,
                                                0})};

        auto& t = ff.getTensor();
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 8; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real());
                EXPECT_NEAR(comRe, 1., 1e-4);
            }
        }
    }

    TEST(FFBtoD2starBLRVarTest, zeroRecoil) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD2starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD2STARPLUS);
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

    TEST(FFBtoD2starBLRVarTest, uninitialized) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "GeV");

        FFBtoD2starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto& pdg = PID::instance();
        const double Mb = pdg.getMass(PID::BZERO);
        const double Mc = pdg.getMass(PID::DSSD2STARPLUS);
        const double Sqq = (Mb - Mc) * (Mb - Mc) + 1e-6;

        ff.resetInitialized();

        ParseOutput po{};
        po.divert();

        static_cast<FormFactorBase&>(ff).evalAtPSPoint({Sqq}, {Mb, Mc});
        auto tmp = po.getStream().str().substr(po.getStream().str().size() - 41);
        EXPECT_STREQ(tmp.c_str(), "Warning, Settings have not been defined!\n");

        po.restore();
    }

    TEST(FFBtoD2starBLRVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBtoD2starBLRVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();

        auto cloned = ff.clone("testLabel");
        ASSERT_NE(cloned, nullptr);
        EXPECT_NE(cloned->group().find("testLabel"), string::npos);

        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pD2starmes{3189.85994318182, 0, 0, -2029.45447278719};
        FourMomentum kNuTau{50.20101629360, -24.8428481565870, -43.029075211928, -7.17451174540};
        FourMomentum pTau{2039.93904052459, 24.8428481565870, 43.029075211928, 2036.62898453259};

        cloned->setSettingsHandler(set);
        cloned->initSettings();
        cloned->calcUnits();
        cloned->eval(Particle(pBmes, 511), {Particle(pD2starmes, -415), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

        ff.eval(Particle(pBmes, 511), {Particle(pD2starmes, -415), Particle(kNuTau, 14), Particle(pTau, -13)}, {});

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
