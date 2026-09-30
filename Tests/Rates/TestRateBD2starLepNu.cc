///
/// @file  TestRateBD2starLepNu.cc
/// @brief Tests for RateBD2starLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateBD2starLepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateBD2starLepNu.hh"
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

    TEST(RateBD2starLepNuTest, evalAtPSPoint) {

        // Create masses (using defs in Pdg.cc)
        // PID& pdg = PID::instance();
        // double Mb = pdg.getMass(PID::BPLUS);
        // double Mc = pdg.getMass(-PID::DSSD2STAR);
        // double Mt = pdg.getMass(PID::ANTITAU);
        double Sqq = 4.41000000000000;


        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BDSSD2STAR}, {0.8, -1.9, 0.1, 1.4, -1.2, 0.9, 0.3, 0.2})};

        Tensor FFvecConj{"FFvecConj",
                         MD::makeVector({8}, {FF_BDSSD2STAR_HC}, {0.8, -1.9, 0.1, 1.4, -1.2, 0.9, 0.3, 0.2})};

        // Evaluate at Q2 point
        RateBD2starLepNu rateQ2;
        rateQ2.setSignatureIndex(0); // B+ -> neutral D** tau nu
        PID& pdg = PID::instance();
        double Mb = pdg.getMass(PID::BPLUS);
        double Mc = pdg.getMass(-PID::DSSD2STAR);
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_BDSSD2STAR});
        rate.dot(FFvecConj, {FF_BDSSD2STAR_HC});

        // Compare to direct evaluation
        Tensor rateEval{"rateEval", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                                   {0.00598444944333202,
                                                    0.000252279052328884,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    -0.00535094640562696,
                                                    0,
                                                    0.00598444944333202,
                                                    0,
                                                    0.0183731395125036,
                                                    0,
                                                    0.000252279052328884,
                                                    0.000103668861880989,
                                                    0,
                                                    -0.000103668861880989,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    0.000252279052328884,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0.000103668861880989,
                                                    0,
                                                    -0.000103668861880989,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    0.000252279052328884,
                                                    0,
                                                    0,
                                                    -0.000252279052328884,
                                                    -0.000103668861880989,
                                                    0,
                                                    0.000103668861880989,
                                                    0,
                                                    0.000252279052328884,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.000103668861880989,
                                                    0,
                                                    0.000103668861880989,
                                                    0,
                                                    0.000252279052328884,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    0,
                                                    -0.00535094640562696,
                                                    -0.000252279052328884,
                                                    0,
                                                    0.000252279052328884,
                                                    0,
                                                    0.00598444944333202,
                                                    0,
                                                    -0.00535094640562696,
                                                    0,
                                                    -0.0241255926778325,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    0.000252279052328884,
                                                    0,
                                                    0.00598444944333202,
                                                    0,
                                                    -0.00535094640562696,
                                                    0,
                                                    0.0183731395125036,
                                                    0.00598444944333202,
                                                    0.000252279052328884,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    -0.00535094640562696,
                                                    0,
                                                    0.00598444944333202,
                                                    0,
                                                    0.0183731395125036,
                                                    0,
                                                    0,
                                                    0,
                                                    0.000252279052328884,
                                                    0,
                                                    -0.000252279052328884,
                                                    0,
                                                    -0.00535094640562696,
                                                    0,
                                                    0.00598444944333202,
                                                    0,
                                                    -0.0241255926778325,
                                                    0.0183731395125036,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.0241255926778325,
                                                    0,
                                                    0.0183731395125036,
                                                    0,
                                                    0.120380943749096,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0.0183731395125036,
                                                    0,
                                                    -0.0241255926778325,
                                                    0,
                                                    0.120380943749096})};

        rateEval *= pow(GFermi, 2.) / (2. * Mb * Mc); // Convert from dGamma/dw to dGamma/dq^2


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

    TEST(RateBD2starLepNuTest, ranges) {

        RateBD2starLepNu rate;
        vector<pair<double, double>> ranges = {
            {pow(1.77686, 2.), pow(5.27932 - 2.461, 2.)},      {pow(1.77686, 2.), pow(5.27963 - 2.465, 2.)},
            {pow(0.1056583745, 2.), pow(5.27932 - 2.461, 2.)}, {pow(0.1056583745, 2.), pow(5.27963 - 2.465, 2.)},
            {pow(5.11e-4, 2.), pow(5.27932 - 2.461, 2.)},      {pow(5.11e-4, 2.), pow(5.27963 - 2.465, 2.)}};
        for (size_t idx = 0; idx < 6; ++idx) {
            rate.setSignatureIndex(idx);
            rate.init();
            rate.calcTensor();
            const auto& range = rate.getEvaluationPoints();
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
