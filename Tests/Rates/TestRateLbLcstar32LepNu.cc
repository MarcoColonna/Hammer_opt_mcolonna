///
/// @file  TestRateLbLcstar32LepNu.cc
/// @brief Tests for RateLbLcstar32LepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateLbLcstar32LepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateLbLcstar32LepNu.hh"
#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/Pdg.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(RateLbLcstar32LepNuTest, evalAtPSPoint) {

        // Create FF
        Tensor FFvec{"FFvec",
                     MD::makeVector({16}, {FF_LBLCSTAR32},
                                    {1.1, 0.8, 0.9, 0.4, 0.3, 0.1, 1.1, 0.5, 0.2, 0.3, 1.2, 0.3, 0.2, 0.1, 0.5, 0.3})};

        Tensor FFvecConj{"FFvecConj", MD::makeVector({16}, {FF_LBLCSTAR32_HC},
                                                     {1.1, 0.8, 0.9, 0.4, 0.3, 0.1, 1.1, 0.5, 0.2, 0.3, 1.2, 0.3, 0.2,
                                                      0.1, 0.5, 0.3})};

        // Evaluate at FF point
        PID& pdg = PID::instance();
        double Mb = pdg.getMass(PID::LAMBDAB);
        double Mc = pdg.getMass(PID::LAMBDACSTAR32MINUS);
        // double Mt = pdg.getMass(PID::ANTITAU);
        double w = 1.16764653309284;
        double Sqq = pow(Mb, 2.) + pow(Mc, 2.) - 2. * Mb * Mc * w;

        RateLbLcstar32LepNu rateQ2;
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_LBLCSTAR32});
        rate.dot(FFvecConj, {FF_LBLCSTAR32_HC});

        // Compare to direct evaluation
        Tensor rateEval{"", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                           {0.00098183115987019,
                                            0.00098183115987019,
                                            0.00046582921517461,
                                            -0.000195786098358956,
                                            -0.000175409991764773,
                                            -0.00010021682645248,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.00098183115987019,
                                            0.00098183115987019,
                                            0.00046582921517461,
                                            -0.000195786098358956,
                                            -0.000175409991764773,
                                            -0.00010021682645248,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.00046582921517461,
                                            0.00046582921517461,
                                            0.00098183115987019,
                                            -0.000175409991764773,
                                            -0.000195786098358956,
                                            -0.00136536270225896,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.000195786098358956,
                                            -0.000195786098358956,
                                            -0.000175409991764773,
                                            0.0000876880959348696,
                                            0.0000807958732624223,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.000175409991764773,
                                            -0.000175409991764773,
                                            -0.000195786098358956,
                                            0.0000807958732624223,
                                            0.0000876880959348696,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.00010021682645248,
                                            -0.00010021682645248,
                                            -0.00136536270225896,
                                            0,
                                            0,
                                            0.0123122817190737,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.00098183115987019,
                                            0.00046582921517461,
                                            -0.000195786098358956,
                                            -0.000175409991764773,
                                            -0.00136536270225896,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.00046582921517461,
                                            0.00098183115987019,
                                            -0.000175409991764773,
                                            -0.000195786098358956,
                                            -0.00010021682645248,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.000195786098358956,
                                            -0.000175409991764773,
                                            0.0000876880959348696,
                                            0.0000807958732624223,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.000175409991764773,
                                            -0.000195786098358956,
                                            0.0000807958732624223,
                                            0.0000876880959348696,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.00136536270225896,
                                            -0.00010021682645248,
                                            0,
                                            0,
                                            0.0123122817190737})};

        rateEval *= GFermi * GFermi;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                double comRe = compareVals(rate.element({idx1, idx2}).real(), rateEval.element({idx1, idx2}).real());
                double comIm = compareVals(rate.element({idx1, idx2}).imag(), rateEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
