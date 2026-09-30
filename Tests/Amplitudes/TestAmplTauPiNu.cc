///
/// @file  TestAmplTauPiNu.cc
/// @brief Tests for AmplTauPiNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplTauPiNu.hh"
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

    TEST(AmplTauPiNuTest, evalTauDec) {

        // Momentum decs
        FourMomentum pBmes{5280., 0., 0., 0.};
        FourMomentum pDmes{2435.0634286113414, 0., 0., -1560.8884974143486};
        FourMomentum kNuTau{485.1308153380305, 477.3411592184665, 26.61679323421851, 82.39461162260386};

        FourMomentum pTau{2359.8057560506286, -477.3411592184665, -26.61679323421851, 1478.4938857917452};
        FourMomentum kNuBarTau{1624.9382437479949, -1114.3651955416435, -168.38866277860308, 1170.580952009679};
        FourMomentum pPion{734.8675123026333, 637.0240363231768, 141.77186954438457, 307.9129337820657};

        // Evaluate at PS point
        AmplTauPiNu ampl;
        ampl.eval(Particle(pTau, -15), {Particle(pPion, 211), Particle(kNuBarTau, -16)},
                  {Particle(pDmes, -411), Particle(kNuTau, 16)});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval", MD::makeVector({2, 2}, {SPIN_NUTAU_REF, SPIN_TAUP},
                                                   {996617.8037695728, -7.440004860781783e6 - 4.770341272969153e6 * 1i,
                                                    -838974.9989246334 + 537929.3615056705 * 1i, 8.837976475927593e6})};
        amplEval *= GFermi * FPion;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                const Tensor& t = ampl.getTensor();
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

    TEST(AmplTauPiNuTest, evalTauDecNoRef) {

        // Momentum decs
        FourMomentum pDmes{1., 0., 0., 0.};
        FourMomentum kNuTau{1., 0., 0., 1.};
        FourMomentum pTau{2359.8057560506286, -477.3411592184665, -26.61679323421851, 1478.4938857917452};
        // FourMomentum pBmes = pDmes + pTau + pDmes;

        FourMomentum kNuBarTau{1624.9382437479949, -1114.3651955416435, -168.38866277860308, 1170.580952009679};
        FourMomentum pPion{734.8675123026333, 637.0240363231768, 141.77186954438457, 307.9129337820657};

        // Evaluate at PS point
        AmplTauPiNu ampl;
        ampl.eval(Particle(pTau, -15), {Particle(pPion, 211), Particle(kNuBarTau, -16)}, {});

        // Compare to direct evaluation
        Tensor amplEval{"amplEval",
                        MD::makeVector({2, 2}, {SPIN_NUTAU_REF, SPIN_TAUP},
                                       {6.1701224545758795e6, 6.359163619983471e6 - 770520.6947048472 * 1i,
                                        6.125322015933577e6 + 742186.8121424909 * 1i, 6.405674369758315e6})};
        amplEval *= GFermi * FPion;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                const Tensor& t = ampl.getTensor();
                double comRe = compareVals(t.element({idx1, idx2}).real(), amplEval.element({idx1, idx2}).real());
                double comIm = compareVals(t.element({idx1, idx2}).imag(), amplEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
