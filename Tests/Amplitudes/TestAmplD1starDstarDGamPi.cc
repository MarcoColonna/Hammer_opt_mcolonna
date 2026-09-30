///
/// @file  TestAmplD1starDstarDGamPi.cc
/// @brief Tests for AmplD1starDstarDGamPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplD1starDstarDGamPi.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
#include "gtest/gtest.h"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(AmplD1starDstarDGamPiTest, evalTauDec) {

        FourMomentum pDssmes{2430., 0, 0, 0};
        FourMomentum pDsmes{2042.5462962963, 226.188257790144, 261.035164712006, -112.226945674169};
        FourMomentum pPmes{387.453703703704, -226.188257790144, -261.035164712006, 112.226945674169};
        FourMomentum pDmes{1923.84498740737, 253.18202510169, 358.06998615436, -175.04080310966};
        FourMomentum kGmes{118.70130888893, -26.993767311543, -97.034821442352, 62.813857435493};
        FourMomentum kNuTau{103.44530435713, 31.76995444103, 97.77786575997, -11.44945220433};
        FourMomentum pTau{4212.0896750667, -31.76995444103, -97.77786575997, 4209.5254589759};

        // Evaluate at PS point
        AmplD1starDstarDGamPi ampl;
        ampl.eval(Particle(pDssmes, 20423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(kGmes, 22)},
                  {Particle(kNuTau, 16), Particle(pTau, -15)});

        // Compare to direct evaluation
        Tensor amplEval{
            "amplEval",
            MD::makeVector({1, 3, 2}, {FF_DSSD1STAR, SPIN_DSSD1STAR, SPIN_GAMMA},
                           {3.6472794760505e6 + 2.0738197821235e6 * 1i, -1.13715800640331e7 + 5.4445378651272e6 * 1i,
                            -9.1149311284137e6 - 4.76584609262e6 * 1i, -9.1149311284137e6 + 4.76584609262e6 * 1i,
                            -1.13715800640331e7 - 5.4445378651272e6 * 1i, 3.6472794760505e6 - 2.0738197821235e6 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    double comRe =
                        compareVals(t.element({idx1, idx2, idx3}).real(), amplEval.element({idx1, idx2, idx3}).real());
                    double comIm =
                        compareVals(t.element({idx1, idx2, idx3}).imag(), amplEval.element({idx1, idx2, idx3}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }
        }
    }

    TEST(AmplD1starDstarDGamPiTest, evalTauDecNoRef) {

        // Momentum decs

        FourMomentum pDssmes{2430., 0, 0, 0};
        FourMomentum pDsmes{2042.5462962963, 226.188257790144, 261.035164712006, -112.226945674169};
        FourMomentum pPmes{387.453703703704, -226.188257790144, -261.035164712006, 112.226945674169};
        FourMomentum pDmes{1923.84498740737, 253.18202510169, 358.06998615436, -175.04080310966};
        FourMomentum kGmes{118.70130888893, -26.993767311543, -97.034821442352, 62.813857435493};
        //        FourMomentum kNuTau{1., 0.707106781186548, 0, 0.70710678118655};
        //        FourMomentum pTau{1., -0.70710678118655, 0, 0.70710678118655};

        // Evaluate at PS point
        AmplD1starDstarDGamPi ampl;
        ampl.eval(Particle(pDssmes, 20423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(kGmes, 22)}, {});

        // Compare to direct evaluation
        Tensor amplEval{
            "amplEval",
            MD::makeVector({1, 3, 2}, {FF_DSSD1STAR, SPIN_DSSD1STAR, SPIN_GAMMA},
                           {3.0993911587449e6 - 2.8279233565003e6 * 1i, 1.6640517221633e6 + 1.24974700473129e7 * 1i,
                            -9.1149311284137e6 - 4.76584609262e6 * 1i, -9.1149311284137e6 + 4.76584609262e6 * 1i,
                            1.6640517221633e6 - 1.24974700473129e7 * 1i, 3.0993911587449e6 + 2.8279233565003e6 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    double comRe =
                        compareVals(t.element({idx1, idx2, idx3}).real(), amplEval.element({idx1, idx2, idx3}).real());
                    double comIm =
                        compareVals(t.element({idx1, idx2, idx3}).imag(), amplEval.element({idx1, idx2, idx3}).imag());
                    EXPECT_NEAR(comRe, 1., 1e-4);
                    EXPECT_NEAR(comIm, 1., 1e-4);
                }
            }
        }
    }

} // namespace Hammer
