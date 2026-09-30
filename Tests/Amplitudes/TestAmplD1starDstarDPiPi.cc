///
/// @file  TestAmplD1starDstarDPiPi.cc
/// @brief Tests for AmplD1starDstarDPiPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplD1starDstarDPiPi.hh"
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

    TEST(AmplD1starDstarDPiPiTest, evalTauDec) {

        FourMomentum pDssmes{2430., 0, 0, 0};
        FourMomentum pDsmes{2042.5462962963, 226.188257790144, 261.035164712006, -112.226945674169};
        FourMomentum pPmes{387.453703703704, -226.188257790144, -261.035164712006, 112.226945674169};
        FourMomentum pDmes{1902.2671378386, 222.07352102032, 275.03024776467, -124.245634839982};
        FourMomentum pPDmes{140.2791584577, 4.114736769829, -13.995083052665, 12.0186891658129};
        FourMomentum kNuTau{103.44530435713, 31.76995444103, 97.77786575997, -11.44945220433};
        FourMomentum pTau{4212.0896750667, -31.76995444103, -97.77786575997, 4209.5254589759};


        // Evaluate at PS point
        AmplD1starDstarDPiPi ampl;
        ampl.eval(Particle(pDssmes, 20423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(pPDmes, 111)},
                  {Particle(kNuTau, 16), Particle(pTau, -15)});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({1, 3}, {FF_DSSD1STAR, SPIN_DSSD1STAR},
                                                   {-2.7254286501479e7 - 960765.762774 * 1i, 2.2303725367153e7,
                                                    -2.7254286501479e7 + 960765.762774 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

    TEST(AmplD1starDstarDPiPiTest, evalTauDecNoRef) {

        // Momentum decs

        FourMomentum pDssmes{2430., 0, 0, 0};
        FourMomentum pDsmes{2042.5462962963, 226.188257790144, 261.035164712006, -112.226945674169};
        FourMomentum pPmes{387.453703703704, -226.188257790144, -261.035164712006, 112.226945674169};
        FourMomentum pDmes{1902.2671378386, 222.07352102032, 275.03024776467, -124.245634839982};
        FourMomentum pPDmes{140.2791584577, 4.114736769829, -13.995083052665, 12.0186891658129};
        //        FourMomentum kNuTau{1., 0.707106781186548, 0, 0.70710678118655};
        //        FourMomentum pTau{1., -0.70710678118655, 0, 0.70710678118655};


        // Evaluate at PS point
        AmplD1starDstarDPiPi ampl;
        ampl.eval(Particle(pDssmes, 20423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(pPDmes, 111)}, {});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({1, 3}, {FF_DSSD1STAR, SPIN_DSSD1STAR},
                                                   {-9.33578023784e6 + 2.5623473825896e7 * 1i, 2.2303725367153e7,
                                                    -9.33578023784e6 - 2.5623473825896e7 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
