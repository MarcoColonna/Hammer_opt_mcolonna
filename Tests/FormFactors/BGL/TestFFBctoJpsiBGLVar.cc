///
/// @file  TestFFBctoJpsiBGLVar.cc
/// @brief Tests for FFBctoJpsiBGLVar
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#define FORMFACTORBASE_FRIENDS                 \
    FRIEND_TEST(FFBctoJpsiBGLVarTest, eval);   \
    FRIEND_TEST(FFBctoJpsiBGLVarTest, clone);  \
    FRIEND_TEST(FFBctoJpsiBGLVarTest, addRefs)

#include "Hammer/FormFactors/BGL/FFBctoJpsiBGLVar.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/SettingsConsumer.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(FFBctoJpsiBGLVarTest, eval) {

        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");
        // Momentum decs
        FourMomentum pBcmes{6275., 0., 0., 0.};
        FourMomentum pJpsimes{3821.5809561753, -857.35310207508, 0., -2069.83348677229};
        FourMomentum kNuTau{131.0024243058, -0.827618332, 98.867858366606, -85.9424039919};
        FourMomentum pTau{2322.4166195189, 858.1807204071, -98.867858366606, 2155.7758907642};

        // Evaluate the FF class
        FFBctoJpsiBGLVar ff;
        //_mFFErrNames = ;
        ff.setSettingsHandler(set);
        ff.initSettings();
        ff.calcUnits();
        ff.eval(Particle(pBcmes, 541), {Particle(pJpsimes, 443), Particle(kNuTau, 12), Particle(pTau, -11)}, {});


        // Compare to direct evaluation
        Tensor ffEval{"ffEval", MD::makeVector({8, 16}, {FF_BCJPSI, FF_BCJPSI_VAR},
                                               {0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                4427.03290506073,
                                                -0.393559090703304,
                                                5.2200664925078,
                                                107.397456924871,
                                                1.31325024502913,
                                                -1.00446946868195,
                                                18.2164992854631,
                                                -120.181097801326,
                                                -1.58924967192751,
                                                2.45288559054112,
                                                -0.570342580089904,
                                                120.975602889146,
                                                1.48533161940981,
                                                8.09786324021549,
                                                -195.838127251623,
                                                0.,
                                                0.0000826899544190131,
                                                -1.06650890118958e-7,
                                                9.70984284586356e-7,
                                                3.44807109596066e-9,
                                                1.17621647364217e-7,
                                                -5.84871082482626e-7,
                                                1.84426858205508e-6,
                                                1.92010490295894e-7,
                                                1.6471584872917e-6,
                                                3.63453097610353e-7,
                                                7.61893411650106e-7,
                                                -2.92774398341249e-10,
                                                -1.19709315552788e-6,
                                                5.45964640343237e-6,
                                                3.0343254229567e-8,
                                                0.,
                                                0.00025086451858624,
                                                -8.01760378796252e-6,
                                                1.51425780677525e-6,
                                                0.0000114213836970665,
                                                3.25932822130267e-6,
                                                0.0000149006926937093,
                                                -9.28181995711475e-6,
                                                -0.0000324332423499511,
                                                -0.0000191685803930914,
                                                5.32782053353738e-6,
                                                -0.0000214358649775044,
                                                0.0000242410838470572,
                                                -0.0000224848024015673,
                                                0.000133700162812169,
                                                0.00013651808761851,
                                                0.,
                                                -0.000048116520554789,
                                                1.00557807714362e-6,
                                                3.70873623940592e-7,
                                                -3.97694930703369e-6,
                                                -8.24183667903914e-8,
                                                1.29917669281939e-6,
                                                1.74303062991995e-6,
                                                5.32751136966085e-6,
                                                -6.02852146093275e-7,
                                                1.59406814125734e-6,
                                                -1.59485532724486e-7,
                                                -4.86344269427056e-6,
                                                4.36565881255317e-6,
                                                -9.13305624775198e-7,
                                                2.01589516006069e-6,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.,
                                                0.})};


        auto& t = ff.getTensor();
        EXPECT_TRUE(t.isSameLabelShape(ffEval));
        for (IndexType idx1 = 1; idx1 < 5; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 15; ++idx2) {
                EXPECT_NEAR(compareVals(t.element({idx1, idx2}).real(), ffEval.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2}).imag(), ffEval.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBctoJpsiBGLVarTest, clone) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBctoJpsiBGLVar ff;
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
        for (IndexType idx1 = 0; idx1 < 8; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 16; ++idx2) {
                EXPECT_NEAR(compareVals(tOrig.element({idx1, idx2}).real(), tCloned.element({idx1, idx2}).real()), 1.,
                            1.e-4);
                EXPECT_DOUBLE_EQ(tOrig.element({idx1, idx2}).imag(), tCloned.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(FFBctoJpsiBGLVarTest, addRefs) {
        SettingsHandler set{};
        set.addSetting<string>("Hammer", "Units", "MeV");

        FFBctoJpsiBGLVar ff;
        ff.setSettingsHandler(set);
        ff.initSettings();

        EXPECT_FALSE(set.checkReference("Cohen:2019zev"));

        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Cohen:2019zev"));

        // Second call must be idempotent.
        static_cast<SettingsConsumer*>(&ff)->addRefs();
        EXPECT_TRUE(set.checkReference("Cohen:2019zev"));

        const string tmpfile = testing::TempDir() + "hammer_bgl_FFBctoJpsiBGLVar_refs_test.bib";
        set.saveReferences(tmpfile);

        ifstream reffile(tmpfile);
        ASSERT_TRUE(reffile.is_open());
        string content((istreambuf_iterator<char>(reffile)), istreambuf_iterator<char>());
        reffile.close();
        std::filesystem::remove(tmpfile);

        EXPECT_NE(content.find("Cohen:2019zev"), string::npos);
    }


} // namespace Hammer
