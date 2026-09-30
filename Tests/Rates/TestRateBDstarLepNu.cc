///
/// @file  TestRateBDstarLepNu.cc
/// @brief Tests for RateBDstarLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateBDstarLepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateBDstarLepNu.hh"
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

    TEST(RateBDstarLepNuTest, evalAtPSPoint) {

        // Create FF and masses (using defs Pdg.cc 5.27932,2.00685,1.77686)
        PID& pdg = PID::instance();
        double Mb = pdg.getMass(PID::BPLUS);
        double Mds = pdg.getMass(-PID::DSTAR);
        double Mt = pdg.getMass(PID::ANTITAU);
        double Sqq = 5.018097451146302;
        double Ew = (Mb * Mb - Mds * Mds + Sqq) / (2 * Mb);
        double Pw = sqrt(Ew * Ew - Sqq);
        // Form Factors Set Ff with model, equivalent to h_A1. Use HQET LO for others.
        double Ff = 1.;
        // Axial and vector FFs
        double MbdsSqq = pow(Mb + Mds, 2.) - Sqq;
        double Fg = Ff / MbdsSqq;
        double Fp = -Ff / MbdsSqq;
        double Fm = Ff / MbdsSqq;
        // NP form factors
        double Fs = -2 * Mds * Ff / MbdsSqq;
        double Fpt = -(Mb + Mds) * Ff / MbdsSqq;
        double Fmt = (Mb - Mds) * Ff / MbdsSqq;
        double Fzt = 0 * Ff;


        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BDSTAR}, {Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt})};

        Tensor FFvecConj{"FFvecConj", MD::makeVector({8}, {FF_BDSTAR_HC}, {Fs, Ff, Fg, Fm, Fp, Fzt, Fmt, Fpt})};

        // Evaluate at FF point
        RateBDstarLepNu rateQ2;
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_BDSTAR});
        rate.dot(FFvecConj, {FF_BDSTAR_HC});

        // Compare to direct evaluation
        Tensor rateEval{"rateEval", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                                   {2.3093000935553447,
                                                    0.1527191272457926,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    -2.101790940990059,
                                                    0.,
                                                    2.3093000935553447,
                                                    0.,
                                                    5.730747726559924,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.05919431698491795,
                                                    0.,
                                                    -0.05919431698491795,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.05919431698491795,
                                                    0.,
                                                    -0.05919431698491795,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.,
                                                    0.,
                                                    -0.15271912724579256,
                                                    -0.05919431698491795,
                                                    0.,
                                                    0.05919431698491795,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -0.05919431698491795,
                                                    0.,
                                                    0.05919431698491795,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    0.,
                                                    -2.101790940990059,
                                                    -0.15271912724579256,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.,
                                                    2.3093000935553447,
                                                    0.,
                                                    -2.101790940990059,
                                                    0.,
                                                    -8.174253762492604,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.,
                                                    2.3093000935553447,
                                                    0.,
                                                    -2.101790940990059,
                                                    0.,
                                                    5.730747726559924,
                                                    2.3093000935553447,
                                                    0.1527191272457926,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    -2.101790940990059,
                                                    0.,
                                                    2.3093000935553447,
                                                    0.,
                                                    5.730747726559924,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.1527191272457926,
                                                    0.,
                                                    -0.15271912724579256,
                                                    0.,
                                                    -2.101790940990059,
                                                    0.,
                                                    2.3093000935553447,
                                                    0.,
                                                    -8.174253762492604,
                                                    5.730747726559924,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    -8.174253762492604,
                                                    0.,
                                                    5.730747726559924,
                                                    0.,
                                                    46.81160147572157,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    0.,
                                                    5.730747726559924,
                                                    0.,
                                                    -8.174253762492604,
                                                    0.,
                                                    46.81160147572157})};
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

    TEST(RateBDstarLepNuTest, ranges) {

        RateBDstarLepNu rate;
        vector<pair<double, double>> ranges = {
            {pow(1.77686, 2.), pow(5.27932 - 2.00685, 2.)},      {pow(1.77686, 2.), pow(5.27963 - 2.01026, 2.)},
            {pow(0.1056583745, 2.), pow(5.27932 - 2.00685, 2.)}, {pow(0.1056583745, 2.), pow(5.27963 - 2.01026, 2.)},
            {pow(5.11e-4, 2.), pow(5.27932 - 2.00685, 2.)},      {pow(5.11e-4, 2.), pow(5.27963 - 2.01026, 2.)}};
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
