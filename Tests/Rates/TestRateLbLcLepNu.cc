///
/// @file  TestRateLbLcLepNu.cc
/// @brief Tests for RateLbLcLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateLbLcLepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateLbLcLepNu.hh"
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

    TEST(RateLbLcLepNuTest, evalAtPSPoint) {

        // Create FF
        Tensor FFvec{"FFvec",
                     MD::makeVector({12}, {FF_LBLC}, {0.8, 0.9, 1.0, 0.3, 0.3, 1.1, 0.4, 0.4, 0.9, 0.2, 0.2, 0.2})};

        Tensor FFvecConj{"FFvecConj", MD::makeVector({12}, {FF_LBLC_HC},
                                                     {0.8, 0.9, 1.0, 0.3, 0.3, 1.1, 0.4, 0.4, 0.9, 0.2, 0.2, 0.2})};

        // Evaluate at FF point
        PID& pdg = PID::instance();
        double Mb = pdg.getMass(PID::LAMBDAB);
        double Mc = pdg.getMass(PID::LAMBDACMINUS);
        double Mt = pdg.getMass(PID::ANTITAU);
        double w = 1.1;
        double Sqq = pow(Mb, 2.) + pow(Mc, 2.) - 2. * Mb * Mc * w;

        RateLbLcLepNu rateQ2;
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_LBLC});
        rate.dot(FFvecConj, {FF_LBLC_HC});

        // Compare to direct evaluation
        Tensor rateEval{"", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                           {10.717384630110057,
                                            -1.7878982028885253,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            -1.9699061510079816,
                                            0.,
                                            10.717384630110057,
                                            0.,
                                            14.724002552109953,
                                            0.,
                                            -1.7878982028885253,
                                            1.4250000000000003,
                                            0.,
                                            1.2630000000000001,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            -1.7878982028885253,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            1.4250000000000003,
                                            0.,
                                            1.2630000000000001,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            -1.7878982028885253,
                                            0.,
                                            0.,
                                            -2.0085556418901187,
                                            1.2630000000000001,
                                            0.,
                                            1.4250000000000003,
                                            0.,
                                            -1.7878982028885253,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            1.2630000000000001,
                                            0.,
                                            1.4250000000000003,
                                            0.,
                                            -1.7878982028885253,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            0.,
                                            -1.9699061510079816,
                                            -2.0085556418901187,
                                            0.,
                                            -1.7878982028885253,
                                            0.,
                                            10.717384630110057,
                                            0.,
                                            -1.9699061510079816,
                                            0.,
                                            -17.654088315154304,
                                            0.,
                                            0.,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            -1.7878982028885253,
                                            0.,
                                            10.717384630110057,
                                            0.,
                                            -1.9699061510079816,
                                            0.,
                                            14.724002552109955,
                                            10.717384630110057,
                                            -1.7878982028885253,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            -1.9699061510079816,
                                            0.,
                                            10.717384630110057,
                                            0.,
                                            14.724002552109953,
                                            0.,
                                            0.,
                                            0.,
                                            -1.7878982028885253,
                                            0.,
                                            -2.0085556418901187,
                                            0.,
                                            -1.9699061510079816,
                                            0.,
                                            10.717384630110057,
                                            0.,
                                            -17.654088315154304,
                                            14.724002552109953,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            -17.654088315154304,
                                            0.,
                                            14.724002552109953,
                                            0.,
                                            59.902621454630584,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            0.,
                                            14.724002552109955,
                                            0.,
                                            -17.654088315154304,
                                            0.,
                                            59.902621454630584})};
        double RateNorm = (pow(GFermi, 2.) * pow(Mb, 5.) * pow(Mc / Mb, 3.) *
                           pow(Sqq / pow(Mb, 2.) - pow(Mt / Mb, 2.), 2.) * sqrt(-1 + pow(w, 2.))) /
                          (16. * Sqq / pow(Mb, 2.) * pi3) / (2. * Mb * Mc);

        rateEval *= RateNorm;

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
