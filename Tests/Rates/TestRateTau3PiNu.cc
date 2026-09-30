///
/// @file  TestRateTau3PiNu.cc
/// @brief Tests for RateTau3PiNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateTau3PiNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateTau3PiNu.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {


    namespace MD = MultiDimensional;

    TEST(RateTau3PiNuTest, evalAtPSPoint) {

        // Using PDG values Mt = 1.77686 GeV, Mp = 0.139571 GeV
        double Sqp = 1.181003641;
        double s1 = 0.2490221257;
        double s2 = 0.1498536045;

        // Evaluate at FF point
        RateTau3PiNu rateP2s1s2;
        Tensor rate = static_cast<RateBase*>(&rateP2s1s2)->evalAtPSPoint({Sqp, s1, s2});

        // Compare to direct evaluation
        Tensor rateEval{"rateEval", MD::makeVector({3, 3}, {FF_TAU3PI, FF_TAU3PI_HC},
                                                   {5.348170407e-8, -6.348410606e-8, 0., -6.348410606e-8,
                                                    9.966642831e-8, 0., 0., 0., 7.423868365e-7})};
        double RateNorm = pow(GFermi, 2.);
        rateEval *= RateNorm;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 3; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                double comRe = compareVals(rate.element({idx1, idx2}).real(), rateEval.element({idx1, idx2}).real());
                double comIm = compareVals(rate.element({idx1, idx2}).imag(), rateEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

    TEST(RateTau3PiNuTest, boundaryFunctions) {
        RateTau3PiNu rate;
        rate.setSignatureIndex(0);
        const IntegrationBoundaries& bounds = rate.getIntegrationBoundaries();

        const double Sqp = 1.181003641;
        const double s1 = 0.2490221257;

        auto s1bounds = bounds[1]({Sqp});
        EXPECT_NEAR(s1bounds.first, 0.07791982, 1e-7);  // 4 * Mp^2
        EXPECT_NEAR(s1bounds.second, 0.89712969, 1e-7); // (sqrt(Sqp) - Mp)^2

        auto s2bounds = bounds[2]({Sqp, s1});
        EXPECT_NEAR(s2bounds.first, 0.12145103, 1e-7);
        EXPECT_NEAR(s2bounds.second, 0.86897035, 1e-7);
    }

} // namespace Hammer
