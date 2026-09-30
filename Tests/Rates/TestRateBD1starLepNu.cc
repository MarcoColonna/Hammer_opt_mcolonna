///
/// @file  TestRateBD1starLepNu.cc
/// @brief Tests for RateBD1starLepNu
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define RATEBASE_FRIENDS FRIEND_TEST(RateBD1starLepNuTest, evalAtPSPoint)

#include "Hammer/Rates/RateBD1starLepNu.hh"
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

    TEST(RateBD1starLepNuTest, evalAtPSPoint) {

        // Create masses (using defs in Pdg.cc)
        // PID& pdg = PID::instance();
        // double Mb = pdg.getMass(PID::BPLUS);
        // double Mc = pdg.getMass(-PID::DSSD1STAR);
        // double Mt = pdg.getMass(PID::ANTITAU);
        double Sqq = 4.41000000000000;


        Tensor FFvec{"FFvec", MD::makeVector({8}, {FF_BDSSD1STAR}, {0.8, 0.5, 0.1, -1.1, 1.3, -0.9, 0.8, 0.2})};

        Tensor FFvecConj{"FFvecConj",
                         MD::makeVector({8}, {FF_BDSSD1STAR_HC}, {0.8, 0.5, 0.1, -1.1, 1.3, -0.9, 0.8, 0.2})};

        // Evaluate at Q2 point
        RateBD1starLepNu rateQ2;
        rateQ2.setSignatureIndex(0); // B+ -> neutral D** tau nu
        PID& pdg = PID::instance();
        double Mb = pdg.getMass(PID::BPLUS);
        double Mc = pdg.getMass(-PID::DSSD1STAR);
        Tensor rate = static_cast<RateBase*>(&rateQ2)->evalAtPSPoint({Sqq});
        rate.dot(FFvec, {FF_BDSSD1STAR});
        rate.dot(FFvecConj, {FF_BDSSD1STAR_HC});

        // Compare to direct evaluation
        Tensor rateEval{"rateEval", MD::makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                                                   {0.00420581638459869,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00126076368898061,
                                                    0,
                                                    0.00420581638459869,
                                                    0,
                                                    0.00474973851639416,
                                                    0,
                                                    -0.00031448741002627,
                                                    0.000571679203575662,
                                                    0,
                                                    0.000571679203575662,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0.000571679203575662,
                                                    0,
                                                    0.000571679203575662,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    0,
                                                    -0.00031448741002627,
                                                    0.000571679203575662,
                                                    0,
                                                    0.000571679203575662,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0.000571679203575662,
                                                    0,
                                                    0.000571679203575662,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    0,
                                                    -0.00126076368898061,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    0.00420581638459869,
                                                    0,
                                                    -0.00126076368898061,
                                                    0,
                                                    -0.0105154027497912,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    0.00420581638459869,
                                                    0,
                                                    -0.00126076368898061,
                                                    0,
                                                    0.00474973851639416,
                                                    0.00420581638459869,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00126076368898061,
                                                    0,
                                                    0.00420581638459869,
                                                    0,
                                                    0.00474973851639416,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00031448741002627,
                                                    0,
                                                    -0.00126076368898061,
                                                    0,
                                                    0.00420581638459869,
                                                    0,
                                                    -0.0105154027497912,
                                                    0.00474973851639416,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    -0.0105154027497912,
                                                    0,
                                                    0.00474973851639416,
                                                    0,
                                                    0.0320791449866136,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0.00474973851639416,
                                                    0,
                                                    -0.0105154027497912,
                                                    0,
                                                    0.0320791449866136})};

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

    TEST(RateBD1starLepNuTest, ranges) {

        RateBD1starLepNu rate;
        vector<pair<double, double>> ranges = {
            {pow(1.77686, 2.), pow(5.27932 - 2.427, 2.)},      {pow(1.77686, 2.), pow(5.27963 - 2.427, 2.)},
            {pow(0.1056583745, 2.), pow(5.27932 - 2.427, 2.)}, {pow(0.1056583745, 2.), pow(5.27963 - 2.427, 2.)},
            {pow(5.11e-4, 2.), pow(5.27932 - 2.427, 2.)},      {pow(5.11e-4, 2.), pow(5.27963 - 2.427, 2.)}};
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
