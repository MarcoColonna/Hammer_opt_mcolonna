///
/// @file  TestRateD2starDstarPi.cc
/// @brief Tests for RateD2starDstarPi
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateD2starDstarPiTest, evalAtPSPoint)

#include "Hammer/Rates/RateD2starDstarPi.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(RateD2starDstarPiTest, evalAtPSPoint) {

        RateD2starDstarPi rate;
        rate.setSignatureIndex(0); // D2*(2460)0 -> D*- pi+
        Tensor result = static_cast<RateBase*>(&rate)->evalAtPSPoint({});

        // Mdss=2.461 GeV, Mds=2.01026 GeV, Mp=0.13957061 GeV, FPion=0.093 GeV
        // Pp^5 / (40 * pi * Mdss^2 * FPion^2) = 1.3557443282e-03
        Tensor rateEval{"", MD::makeVector({1, 1}, {FF_DSSD2STAR, FF_DSSD2STAR_HC}, {1.3557443282e-03})};

        for (IndexType idx1 = 0; idx1 < 1; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 1; ++idx2) {
                double comRe = compareVals(result.element({idx1, idx2}).real(), rateEval.element({idx1, idx2}).real());
                double comIm = compareVals(result.element({idx1, idx2}).imag(), rateEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
