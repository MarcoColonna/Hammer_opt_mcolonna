///
/// @file  TestAmplTau3PiNu.cc
/// @brief Tests for AmplTau3PiNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplTau3PiNu.hh"
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

    TEST(AmplTau3PiNuTest, evalTauDec) {

        // Momentum decs
        //        FourMomentum pBmes{7.27993,	    4.69979,    1.60853,	0.666025};
        //        FourMomentum pDmes{5.42945, 	4.71853,   	1.60413,   	1.07146};
        //        FourMomentum kNuTau{0.0251286,	-0.018733,	0.00440374,	0.0161594};
        //
        //        FourMomentum pTau{1.82536,	    0.0,	    0.0,          -0.421598};
        //        FourMomentum kNuBarTau{0.421598,0.0,        0.0,	        -0.421598};
        //
        //        FourMomentum pPiPlus1{0.554558,	-0.525713,	0.102531,	0.03643};
        //        FourMomentum pPiPlus2{0.360269,	0.169702,	-0.0979545,	0.268475};
        //        FourMomentum pPiMinus{0.48893,	0.356011,	-0.00457666,	-0.304905};

        FourMomentum pBmes{7.27993, 4.69979, 0.666025, -1.60853};
        FourMomentum pDmes{
            5.42945,
            4.71853,
            1.07146,
            -1.60413,
        };
        FourMomentum kNuTau{0.0251286, -0.018733, 0.0161594, -0.00440374};
        FourMomentum pTau{1.82536, 0.0, -0.421598, 0.0};

        FourMomentum kNuBarTau{0.421598, 0.0, -0.421598, 0.0};
        FourMomentum pPiPlus1{0.554558, -0.525713, 0.03643, -0.102531};
        FourMomentum pPiPlus2{0.360269, 0.169702, 0.268475, 0.0979545};
        FourMomentum pPiMinus{0.48893, 0.356011, -0.304905, 0.00457666};

        // Evaluate at PS point
        AmplTau3PiNu ampl;
        ampl.eval(
            Particle(pTau, -15),
            {Particle(pPiPlus1, 211), Particle(pPiPlus2, 211), Particle(pPiMinus, -211), Particle(kNuBarTau, -16)},
            {Particle(pDmes, -411), Particle(kNuTau, 16)});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval",
                        MD::makeVector({2, 2, 3}, {SPIN_NUTAU_REF, SPIN_TAUP, FF_TAU3PI},
                                       {-1.0449 - 0.140778 * 1i, 0.154839 - 0.19279 * 1i, 0.944604,
                                        0.930928 + 0.575143 * 1i, 0.9561 + 0.610998 * 1i, 2.24279 + 1.24279 * 1i,
                                        -0.982199 + 0.383314 * 1i, 0.0419926 - 0.24368 * 1i, 0.826235 - 0.457838 * 1i,
                                        1.09303 + 0.0518621 * 1i, 1.13243 + 0.0710231 * 1i, 2.5641})};

        //        {4.487233821758619e6 - 1.7703273610403342e6 * 1i, -9.763265665153343e6 - 4.728473180668875e6 * 1i,
        //             -3.494330545604125e6 + 3.325504529306864e6 * 1i, 1.0824307593652954e7 + 717063.5861999617 * 1i});

        amplEval *= GFermi;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 3; ++idx3) {
                    const Tensor& t = ampl.getTensor();
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

    TEST(AmplTau3PiNuTest, evalTauDecNoRef) {

        // Momentum decs
        FourMomentum pDmes{1., 0., 0., 0.};
        FourMomentum kNuTau{1., 0., 0., 1.};
        FourMomentum pTau{1.82536, 0.0, -0.421598, 0.0};
        // FourMomentum pBmes = pDmes + pTau + pDmes;

        FourMomentum kNuBarTau{0.421598, 0.0, -0.421598, 0.0};
        FourMomentum pPiPlus1{0.554558, -0.525713, 0.03643, -0.102531};
        FourMomentum pPiPlus2{0.360269, 0.169702, 0.268475, 0.0979545};
        FourMomentum pPiMinus{0.48893, 0.356011, -0.304905, 0.00457666};

        // Evaluate at PS point
        AmplTau3PiNu ampl;
        ampl.eval(
            Particle(pTau, -15),
            {Particle(pPiPlus1, 211), Particle(pPiPlus2, 211), Particle(pPiMinus, -211), Particle(kNuBarTau, -16)}, {});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval",
                        MD::makeVector({2, 2, 3}, {SPIN_NUTAU_REF, SPIN_TAUP, FF_TAU3PI},
                                       {0.541304 + 1.06431 * 1i, 0.579399 + 0.224889 * 1i, 1.69445,
                                        0.419085 - 0.841233 * 1i, 0.964734 - 0.177753 * 1i, 2.14377 + 0. * 1i,
                                        0.541304 + 1.06431 * 1i, 0.579399 + 0.224889 * 1i, 1.69445 + 0. * 1i,
                                        0.419085 - 0.841233 * 1i, 0.964734 - 0.177753 * 1i, 2.14377})};

        amplEval *= GFermi;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 3; ++idx3) {
                    const Tensor& t = ampl.getTensor();
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
