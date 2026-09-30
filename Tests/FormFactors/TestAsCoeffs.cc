///
/// @file  TestAsCoeffs.cc
/// @brief Tests for AsCoeffs
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"
#include "Hammer/FormFactors/AsCoeffs.hh"
#include "Hammer/Math/Utils.hh"
// May need other stuff here

using namespace std;
using namespace complex_literals;

namespace Hammer {

    TEST(AsCoeffs, Cs) {


        double w = 1.274123699;
        double z = 0.2781316348;

        // Evaluate the FF class
        AsCoeffs C;
        vector<double> Cexp = {C.CS(w, z),  C.CP(w, z),  C.CV1(w, z), C.CV2(w, z), C.CV3(w, z), C.CA1(w, z),
                               C.CA2(w, z), C.CA3(w, z), C.CT1(w, z), C.CT2(w, z), C.CT3(w, z)};

        // Compare to direct evaluation
        vector<double> CexpEval = {-0.5881726190, 0.6352874883,  0.9750993329, -0.4168549682,
                                   -0.1870099718, -0.2483607744, -1.122430616, 0.3771921162,
                                   0.5315765071,  -0.8555009276, 0.3679591796};

        for (size_t idx1 = 0; idx1 < 11; ++idx1) {
            double comRe = compareVals(Cexp[idx1], CexpEval[idx1]);
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(AsCoeffs, CsZeroRec) {

        double w = 1.;
        double z = 0.2781316348;

        // Evaluate the FF class
        AsCoeffs C;
        vector<double> Cexp = {C.CS(w, z),  C.CP(w, z),  C.CV1(w, z), C.CV2(w, z), C.CV3(w, z), C.CA1(w, z),
                               C.CA2(w, z), C.CA3(w, z), C.CT1(w, z), C.CT2(w, z), C.CT3(w, z)};

        // Compare to direct evaluation
        vector<double> CexpEval = {-2. / 3.,
                                   2. / 3.,
                                   0.932418972668,
                                   -0.46818714407,
                                   -0.1984795226,
                                   -0.400914360665,
                                   -1.2154359070481,
                                   0.4063130426509,
                                   0.35433640800177,
                                   2. * (-0.46818714407),
                                   -2. * (-0.1984795226)};

        for (size_t idx1 = 0; idx1 < 11; ++idx1) {
            double comRe = compareVals(Cexp[idx1], CexpEval[idx1]);
            EXPECT_NEAR(comRe, 1., 1e-4);
        }

        w = 0.9;

        // Evaluate the FF class
        Cexp = {C.CS(w, z),  C.CP(w, z),  C.CV1(w, z), C.CV2(w, z), C.CV3(w, z), C.CA1(w, z),
                C.CA2(w, z), C.CA3(w, z), C.CT1(w, z), C.CT2(w, z), C.CT3(w, z)};

        for (size_t idx1 = 0; idx1 < 11; ++idx1) {
            double comRe = compareVals(Cexp[idx1], CexpEval[idx1]);
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

    TEST(AsCoeffs, derCs) {


        double w = 1.274123699;
        double z = 0.2781316348;

        // Evaluate the FF class
        AsCoeffs C;
        vector<double> Cexp = {C.derCS(w, z),  C.derCP(w, z),  C.derCV1(w, z), C.derCV2(w, z),
                               C.derCV3(w, z), C.derCA1(w, z), C.derCA2(w, z), C.derCA3(w, z),
                               C.derCT1(w, z), C.derCT2(w, z), C.derCT3(w, z)};

        // Compare to direct evaluation
        vector<double> CexpEval = {0.21962949798381948, -0.14211154654233704, 0.1154628603073411,
                                   0.1672320273464558,  0.038448791861216636, 0.47720390483349784,
                                   0.30710477873202047, -0.09647035180585362, 0.5630620404500565,
                                   0.26572493786592943, -0.0960161066602272};

        for (size_t idx1 = 0; idx1 < 11; ++idx1) {
            double comRe = compareVals(Cexp[idx1], CexpEval[idx1]);
            EXPECT_NEAR(comRe, 1., 1e-4);
        }
    }

} // namespace Hammer
