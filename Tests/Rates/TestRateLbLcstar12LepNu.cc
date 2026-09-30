///
/// @file  TestRateLbLcstar12LepNu.cc
/// @brief Tests for RateLbLcstar12LepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateLbLcstar12LepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateLbLcstar12LepNu.hh"
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

    TEST(RateLbLcstar12LepNuTest, evalAtPSPoint) {

        // Create FF
        Tensor FFvec{"FFvec", MD::makeVector({12}, {FF_LBLCSTAR12},
                                             {1.1, 0.8, 0.9, 1.3, 0.3, 1.1, 0.8, 0.2, 1.2, 0.8, 0.2, 0.1})};

        Tensor FFvecConj{"FFvecConj", MD::makeVector({12}, {FF_LBLCSTAR12_HC},
                                                     {1.1, 0.8, 0.9, 1.3, 0.3, 1.1, 0.8, 0.2, 1.2, 0.8, 0.2, 0.1})};

        // Evaluate at FF point
        PID& pdg = PID::instance();
        double Mb = pdg.getMass(PID::LAMBDAB);
        double Mc = pdg.getMass(PID::LAMBDACSTAR12MINUS);
        // double Mt = pdg.getMass(PID::ANTITAU);
        double w = 1.17741366811652;
        double Sqq = pow(Mb, 2.) + pow(Mc, 2.) - 2. * Mb * Mc * w;

        RateLbLcstar12LepNu rateQ2;
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_LBLCSTAR12});
        rate.dot(FFvecConj, {FF_LBLCSTAR12_HC});

        // Compare to direct evaluation
        Tensor rateEval{"", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                           {0.00434414201927587,
                                            0.00434414201927587,
                                            -0.00294669728792324,
                                            -0.000669722874981501,
                                            0.000604139583188042,
                                            -0.005709908961061,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.00434414201927587,
                                            0.00434414201927587,
                                            -0.00294669728792324,
                                            -0.000669722874981501,
                                            0.000604139583188042,
                                            -0.005709908961061,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.00294669728792324,
                                            -0.00294669728792324,
                                            0.00434414201927587,
                                            0.000604139583188042,
                                            -0.000669722874981501,
                                            -0.0104207625230386,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.000669722874981501,
                                            -0.000669722874981501,
                                            0.000604139583188042,
                                            0.000213844910180711,
                                            -0.000156755287409104,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.000604139583188042,
                                            0.000604139583188042,
                                            -0.000669722874981501,
                                            -0.000156755287409104,
                                            0.000213844910180711,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.005709908961061,
                                            -0.005709908961061,
                                            -0.0104207625230386,
                                            0,
                                            0,
                                            0.149244019230084,
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
                                            0.00434414201927587,
                                            -0.00294669728792324,
                                            -0.000669722874981501,
                                            0.000604139583188042,
                                            -0.0104207625230386,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.00294669728792324,
                                            0.00434414201927587,
                                            0.000604139583188042,
                                            -0.000669722874981501,
                                            -0.005709908961061,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.000669722874981501,
                                            0.000604139583188042,
                                            0.000213844910180711,
                                            -0.000156755287409104,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.000604139583188042,
                                            -0.000669722874981501,
                                            -0.000156755287409104,
                                            0.000213844910180711,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            -0.0104207625230386,
                                            -0.005709908961061,
                                            0,
                                            0,
                                            0.149244019230084})};

        rateEval *= GFermi * GFermi;

        // Check TEs.
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                double comRe = compareVals(rate.element({idx1, idx2}).real(), rateEval.element({idx1, idx2}).real());
                // cout << rate.element({idx1, idx2}).real() << " " << rateEval.element({idx1, idx2}).real() << '\n';
                double comIm = compareVals(rate.element({idx1, idx2}).imag(), rateEval.element({idx1, idx2}).imag());
                EXPECT_NEAR(comRe, 1., 1e-4);
                EXPECT_NEAR(comIm, 1., 1e-4);
            }
        }
    }

} // namespace Hammer
