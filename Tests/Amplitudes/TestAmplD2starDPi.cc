///
/// @file  TestAmplD2starDPi.cc
/// @brief Tests for AmplD2starDPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplD2starDPi.hh"
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

    TEST(AmplD2starDPiTest, evalTauDec) {

        FourMomentum pDssmes{2461., 0, 0, 0};
        FourMomentum pDmes{1933.46627387241, 317.61262700106, 366.54451125733, -157.588618366517};
        FourMomentum pPmes{527.53372612759, -317.61262700106, -366.54451125733, 157.588618366517};
        FourMomentum kNuTau{103.70179046789, 31.76995444103, 97.77786575997, -13.57277821547};
        FourMomentum pTau{4126.6677747495, -31.76995444103, -97.77786575997, 4124.0504461171};


        // Evaluate at PS point
        AmplD2starDPi ampl;
        ampl.eval(Particle(pDssmes, 425), {Particle(pDmes, 411), Particle(pPmes, -211)},
                  {Particle(kNuTau, 16), Particle(pTau, -15)});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval",
                        MD::makeVector({1, 5}, {FF_DSSD2STAR, SPIN_DSSD2STAR},
                                       {81971.238912214 + 84346.410845473 * 1i, -70403.069091315 - 29752.585440155 * 1i,
                                        75756.314078381, -70403.069091315 + 29752.585440155 * 1i,
                                        81971.238912214 - 84346.410845473 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

    TEST(AmplD2starDPiTest, evalTauDecNoRef) {

        // Momentum decs

        FourMomentum pDssmes{2461., 0, 0, 0};
        FourMomentum pDmes{1933.46627387241, 317.61262700106, 366.54451125733, -157.588618366517};
        FourMomentum pPmes{527.53372612759, -317.61262700106, -366.54451125733, 157.588618366517};
        //        FourMomentum kNuTau{1., 0.707106781186548, 0, 0.707106781186548};
        //        FourMomentum pTau{1., -0.70710678118655, 0, 0.707106781186548};

        // Evaluate at PS point
        AmplD2starDPi ampl;
        ampl.eval(Particle(pDssmes, 425), {Particle(pDmes, 411), Particle(pPmes, -211)}, {});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({1, 5}, {FF_DSSD2STAR, SPIN_DSSD2STAR},
                                                   {-16738.548951179 - 116419.165133262 * 1i,
                                                    -50052.135064857 + 57763.243098873 * 1i, 75756.314078381,
                                                    -50052.135064857 - 57763.243098873 * 1i,
                                                    -16738.548951179 + 116419.165133262 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 5; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
