///
/// @file  TestRateBRhoOmegaLepNu.cc
/// @brief Tests for RateBRhoOmegaLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateBRhoOmegaLepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateBRhoOmegaLepNu.hh"
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

    TEST(RateBRhoOmegaLepNuTest, evalAtPSPoint) {

        double Sqq = 1.157;
        
        //Fix sign convention in T2 and T23 definitions, FF indices 6 and 7. Amounts to switching signs in FFvec for purposes of this test.
        const double signfix = -1.;
        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BOMEGA}, {0.7, 1.1, 0.9, 0.6, -0.4, 1.1, signfix*0.3, signfix*0.8})};

        Tensor FFvecConj{"FFvecConj", MD::makeVector({8}, {FF_BOMEGA_HC}, {0.7, 1.1, 0.9, 0.6, -0.4, 1.1, signfix*0.3, signfix*0.8})};

        // Evaluate at FF point
        RateBRhoOmegaLepNu rateQ2;
        rateQ2.setSignatureIndex(7); // B -> omega mu nu
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_BOMEGA});
        rate.dot(FFvecConj, {FF_BOMEGA_HC});

        // Compare to direct evaluation
        Tensor rateEval{"rateEval", MD::makeVector({11, 11}, {WILSON_BUTAUNU, WILSON_BUTAUNU_HC},
                                                   {0.0038835824494465,
                                                    -0.000631089326039775,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    -0.00236864737842745,
                                                    0,
                                                    0.0038835824494465,
                                                    0,
                                                    -0.00278718146659662,
                                                    0,
                                                    -0.000631089326039775,
                                                    0.00343382544172077,
                                                    0,
                                                    -0.00343382544172077,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    -0.000631089326039775,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0.00343382544172077,
                                                    0,
                                                    -0.00343382544172077,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    -0.000631089326039775,
                                                    0,
                                                    0,
                                                    0.000631089326039775,
                                                    -0.00343382544172077,
                                                    0,
                                                    0.00343382544172077,
                                                    0,
                                                    -0.000631089326039775,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.00343382544172077,
                                                    0,
                                                    0.00343382544172077,
                                                    0,
                                                    -0.000631089326039775,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    0,
                                                    -0.00236864737842745,
                                                    0.000631089326039775,
                                                    0,
                                                    -0.000631089326039775,
                                                    0,
                                                    0.0038835824494465,
                                                    0,
                                                    -0.00236864737842745,
                                                    0,
                                                    -0.00222053164832356,
                                                    0,
                                                    0,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    -0.000631089326039775,
                                                    0,
                                                    0.0038835824494465,
                                                    0,
                                                    -0.00236864737842745,
                                                    0,
                                                    -0.00278718146659662,
                                                    0.0038835824494465,
                                                    -0.000631089326039775,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    -0.00236864737842745,
                                                    0,
                                                    0.0038835824494465,
                                                    0,
                                                    -0.00278718146659662,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.000631089326039775,
                                                    0,
                                                    0.000631089326039775,
                                                    0,
                                                    -0.00236864737842745,
                                                    0,
                                                    0.0038835824494465,
                                                    0,
                                                    -0.00222053164832356,
                                                    -0.00278718146659662,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.00222053164832356,
                                                    0,
                                                    -0.00278718146659662,
                                                    0,
                                                    0.211793889236537,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.00278718146659662,
                                                    0,
                                                    -0.00222053164832356,
                                                    0,
                                                    0.211793889236537})};

        rateEval *= pow(GFermi, 2.);


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

    TEST(RateBRhoOmegaLepNuTest, ranges) {

        RateBRhoOmegaLepNu rate;
        vector<pair<double, double>> ranges = {
            {pow(1.77686, 2.), pow(5.27932 - 0.77526, 2.)},      {pow(1.77686, 2.), pow(5.27963 - 0.77511, 2.)},
            {pow(0.1056583745, 2.), pow(5.27932 - 0.77526, 2.)}, {pow(0.1056583745, 2.), pow(5.27963 - 0.77511, 2.)},
            {pow(5.11e-4, 2.), pow(5.27932 - 0.77526, 2.)},      {pow(5.11e-4, 2.), pow(5.27963 - 0.77511, 2.)},
            {pow(1.77686, 2.), pow(5.27932 - 0.78265, 2.)},      {pow(0.1056583745, 2.), pow(5.27932 - 0.78265, 2.)},
            {pow(5.11e-4, 2.), pow(5.27932 - 0.78265, 2.)}};
        for (size_t idx = 0; idx < 9; ++idx) {
            rate.setSignatureIndex(idx);
            rate.init();
            rate.calcTensor();
            const auto& range = rate.getEvaluationPoints();
            //            cout << idx << " " << ranges[idx].first << " " << ranges[idx].second << '\n';
            //            cout << idx << " " << range.first << " " << range.second << '\n';
            EXPECT_NEAR(range[9][0],
                        0.5 * ((ranges[idx].first + ranges[idx].second) +
                               0.9924070086164836 * (ranges[idx].second - ranges[idx].first)),
                        1.e-5);
            EXPECT_NEAR(range[18][0],
                        0.5 * ((ranges[idx].first + ranges[idx].second) -
                               0.9924070086164836 * (ranges[idx].second - ranges[idx].first)),
                        1.e-5);
        }
    }

} // namespace Hammer
