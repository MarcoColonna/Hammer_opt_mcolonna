///
/// @file  TestAmplD1DstarDPiPi.cc
/// @brief Tests for AmplD1DstarDPiPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplD1DstarDPiPi.hh"
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

    TEST(AmplD1DstarDPiPiTest, evalTauDec) {

        FourMomentum pDssmes{2420., 0, 0, 0};
        FourMomentum pDsmes{2040.96590909091, 220.585259460554, 254.568959940215, -109.44692784614};
        FourMomentum pPmes{379.034090909091, -220.585259460554, -254.568959940215, 109.44692784614};
        FourMomentum pDmes{1900.63599375279, 216.85554717276, 269.00838493206, -121.656653107506};
        FourMomentum pPDmes{140.329915338116, 3.729712287798, -14.43942499184, 12.2097252613656};
        FourMomentum kNuTau{103.37185812991, 31.76995444103, 97.77786575997, -10.76568694837};
        FourMomentum pTau{4240.0165716222, -31.76995444103, -97.77786575997, 4237.4692549437};


        // Evaluate at PS point
        AmplD1DstarDPiPi ampl;
        ampl.eval(Particle(pDssmes, 10423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(pPDmes, 111)},
                  {Particle(kNuTau, 16), Particle(pTau, -15)});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 3}, {FF_DSSD1, SPIN_DSSD1},
                                                   {4.102676798425e6 + 144627.2095011 * 1i, -3.357452655284e6,
                                                    4.102676798425e6 - 144627.2095011 * 1i,
                                                    -2.4706253138945e6 - 1.5738211372155e6 * 1i, 795541.1050564,
                                                    -2.4706253138945e6 + 1.5738211372155e6 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

    TEST(AmplD1DstarDPiPiTest, evalTauDecNoRef) {

        // Momentum decs

        FourMomentum pDssmes{2420., 0, 0, 0};
        FourMomentum pDsmes{2040.96590909091, 220.585259460554, 254.568959940215, -109.44692784614};
        FourMomentum pPmes{379.034090909091, -220.585259460554, -254.568959940215, 109.44692784614};
        FourMomentum pDmes{1900.63599375279, 216.85554717276, 269.00838493206, -121.656653107506};
        FourMomentum pPDmes{140.329915338116, 3.729712287798, -14.43942499184, 12.2097252613656};
        // FourMomentum kNuTau{1.,1./sqrt2,0.,1./sqrt2};
        // FourMomentum pTau{1.,-1./sqrt2,0.,1./sqrt2};

        // Evaluate at PS point
        AmplD1DstarDPiPi ampl;
        ampl.eval(Particle(pDssmes, 10423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(pPDmes, 111)}, {});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 3}, {FF_DSSD1, SPIN_DSSD1},
                                                   {1.4053455031707e6 - 3.8571852378102e6 * 1i, -3.357452655284e6,
                                                    1.4053455031707e6 + 3.8571852378102e6 * 1i,
                                                    -2.2602580567581e6 + 1.863366826597e6 * 1i, 795541.1050564,
                                                    -2.2602580567581e6 - 1.863366826597e6 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
