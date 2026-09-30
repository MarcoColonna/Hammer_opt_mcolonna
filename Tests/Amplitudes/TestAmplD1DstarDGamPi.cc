///
/// @file  TestAmplD1DstarDGamPi.cc
/// @brief Tests for AmplD1DstarDGamPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Amplitudes/AmplD1DstarDGamPi.hh"
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

    TEST(AmplD1DstarDGamPiTest, evalTauDec) {

        FourMomentum pDssmes{2420., 0, 0, 0};
        FourMomentum pDsmes{2040.96590909091, 220.585259460554, 254.568959940215, -109.44692784614};
        FourMomentum pPmes{379.034090909091, -220.585259460554, -254.568959940215, 109.44692784614};
        FourMomentum pDmes{1921.79653549978, 247.90540978984, 351.9804474643, -172.42272547057};
        FourMomentum kGmes{119.169373591134, -27.320150329288, -97.411487524087, 62.975797624428};
        FourMomentum kNuTau{103.37185812991, 31.76995444103, 97.77786575997, -10.76568694837};
        FourMomentum pTau{4240.0165716222, -31.76995444103, -97.77786575997, 4237.4692549437};

        // Evaluate at PS point
        AmplD1DstarDGamPi ampl;
        ampl.eval(Particle(pDssmes, 10423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(kGmes, 22)},
                  {Particle(kNuTau, 16), Particle(pTau, -15)});

        // Compare to direct evaluation
        Tensor amplEval{
            "amplEval",
            MD::makeVector({2, 3, 2}, {FF_DSSD1, SPIN_DSSD1, SPIN_GAMMA},
                           {-549036.8967447 - 312178.8678552 * 1i, 1.7118011028324e6 - 819584.0744608 * 1i,
                            1.3721003651207e6 + 717418.3843703 * 1i, 1.3721003651207e6 - 717418.3843703 * 1i,
                            1.7118011028324e6 + 819584.0744608 * 1i, -549036.8967447 + 312178.8678552 * 1i,
                            -3829.27612148 - 528057.55580948 * 1i, 391417.34481735 + 150803.24170416 * 1i,
                            457366.78837358 + 450652.35253775 * 1i, 457366.78837358 - 450652.35253775 * 1i,
                            391417.34481735 - 150803.24170416 * 1i, -3829.27612148 + 528057.55580948 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
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

    TEST(AmplD1DstarDGamPiTest, evalTauDecNoRef) {

        // Momentum decs

        FourMomentum pDssmes{2420., 0, 0, 0};
        FourMomentum pDsmes{2040.96590909091, 220.585259460554, 254.568959940215, -109.44692784614};
        FourMomentum pPmes{379.034090909091, -220.585259460554, -254.568959940215, 109.44692784614};
        FourMomentum pDmes{1921.79653549978, 247.90540978984, 351.9804474643, -172.42272547057};
        FourMomentum kGmes{119.169373591134, -27.320150329288, -97.411487524087, 62.975797624428};
        // FourMomentum kNuTau{1.,1./sqrt2,0.,1./sqrt2};
        // FourMomentum pTau{1.,-1./sqrt2,0.,1./sqrt2};

        // Evaluate at PS point
        AmplD1DstarDGamPi ampl;
        ampl.eval(Particle(pDssmes, 10423),
                  {Particle(pDsmes, 413), Particle(pPmes, -211), Particle(pDmes, 411), Particle(kGmes, 22)}, {});

        // Compare to direct evaluation
        Tensor amplEval{
            "amplEval",
            MD::makeVector({2, 3, 2}, {FF_DSSD1, SPIN_DSSD1, SPIN_GAMMA},
                           {-466561.4781563 + 425696.5428836 * 1i, -250495.1429027 - 1.8812850007774e6 * 1i,
                            1.3721003651207e6 + 717418.3843703 * 1i, 1.3721003651207e6 - 717418.3843703 * 1i,
                            -250495.1429027 + 1.8812850007774e6 * 1i, -466561.4781563 - 425696.5428836 * 1i,
                            -503395.89082919 - 159536.9007452 * 1i, 264377.01714285 - 325659.25188607 * 1i,
                            457366.78837358 + 450652.35253775 * 1i, 457366.78837358 - 450652.35253775 * 1i,
                            264377.01714285 + 325659.25188607 * 1i, -503395.89082919 + 159536.9007452 * 1i})};

        const Tensor& t = ampl.getTensor();
        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
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
