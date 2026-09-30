///
/// @file  TestRateBDLepNu.cc
/// @brief Tests for RateBDLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateBDLepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateBDLepNu.hh"
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

    TEST(RateBDLepNuTest, evalAtPSPoint) {

        // Create FF and masses (using defs Pdg.cc 5.27932,1.86483,1.77686)
        PID& pdg = PID::instance();
        double Mb = pdg.getMass(PID::BPLUS);
        double Md = pdg.getMass(-PID::D0);
        double Mt = pdg.getMass(PID::ANTITAU);
        double Sqq = 9.844667828974893;
        double Ew = (Mb * Mb - Md * Md + Sqq) / (2 * Mb);
        double Pw = sqrt(Ew * Ew - Sqq);
        double Fp = 1.;
        double MbdSqq = pow(Mb + Md, 2) - Sqq;
        double Fz = Fp * MbdSqq / pow(Mb + Md, 2);
        double Fs = Fp * MbdSqq / (Mb + Md);
        double Ft = Fp / (Mb + Md);

        Tensor FFvec{"FFvec", MD::makeVector({4}, {FF_BD}, {Fs, Fz, Fp, Ft})};

        Tensor FFvecConj{"FFvecConj", MD::makeVector({4}, {FF_BD_HC}, {Fs, Fz, Fp, Ft})};

        // Evaluate at FF point
        RateBDLepNu rateQ2;
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_BD});
        rate.dot(FFvecConj, {FF_BD_HC});

        // Compare to direct evaluation
        Tensor rateEval{"rateEval", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                                   {9.249911846773028,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    -3.7759201223246373,
                                                    0.,
                                                    -10.245171925817354,
                                                    16.6242194008283,
                                                    0.,
                                                    16.6242194008283,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    16.6242194008283,
                                                    0.,
                                                    16.6242194008283,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    0.,
                                                    -10.245171925817354,
                                                    16.6242194008283,
                                                    0.,
                                                    16.6242194008283,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    16.6242194008283,
                                                    0.,
                                                    16.6242194008283,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    0.,
                                                    9.249911846773028,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    -3.7759201223246373,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    -3.7759201223246373,
                                                    9.249911846773028,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    -3.7759201223246373,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    -10.245171925817354,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    9.249911846773028,
                                                    0.,
                                                    -3.7759201223246373,
                                                    -3.7759201223246373,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -3.7759201223246373,
                                                    0.,
                                                    -3.7759201223246373,
                                                    0.,
                                                    6.40878577416312,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -3.7759201223246373,
                                                    0.,
                                                    -3.7759201223246373,
                                                    0.,
                                                    6.40878577416312})};
        double RateNorm =
            (pow(GFermi, 2.) * Pw * pow(-pow(Mt, 2.) + Sqq, 2.)) / (32. * pow(Mb, 2.) * pow(pi, 3.) * Sqq);
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

    TEST(RateBDLepNuTest, ranges) {

        RateBDLepNu rate;
        vector<pair<double, double>> ranges = {
            {pow(1.77686, 2.), pow(5.27932 - 1.86483, 2.)},      {pow(1.77686, 2.), pow(5.27963 - 1.86959, 2.)},
            {pow(0.1056583745, 2.), pow(5.27932 - 1.86483, 2.)}, {pow(0.1056583745, 2.), pow(5.27963 - 1.86959, 2.)},
            {pow(5.11e-4, 2.), pow(5.27932 - 1.86483, 2.)},      {pow(5.11e-4, 2.), pow(5.27963 - 1.86959, 2.)}};
        for (size_t idx = 0; idx < 6; ++idx) {
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
