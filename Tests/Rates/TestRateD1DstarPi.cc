///
/// @file  TestRateD1DstarPi.cc
/// @brief Tests for RateD1DstarPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateD1DstarPiTest, evalAtPSPoint)

#include "Hammer/Rates/RateD1DstarPi.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(RateD1DstarPiTest, evalAtPSPoint) {

        // Evaluate
        RateD1DstarPi rateQ2;
        rateQ2.setSignatureIndex(0); // D1^* -> D*- pi^+
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({});

        // Compare to direct evaluation
        Tensor rateEval{
            "", MD::makeVector({2, 2}, {FF_DSSD1, FF_DSSD1_HC}, {0.00648562567127407, 0, 0, 0.00144125014917202})};

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 2; ++idx2) {
                double comRe = compareVals(rate.element({idx1, idx2}).real(), rateEval.element({idx1, idx2}).real());
                double comIm = compareVals(rate.element({idx1, idx2}).imag(), rateEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }


} // namespace Hammer
