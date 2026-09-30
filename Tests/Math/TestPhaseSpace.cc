///
/// @file  TestPhaseSpace.cc
/// @brief Tests for PhaseSpace
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/Constants.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/PhaseSpace.hh"
#include "gtest/gtest.h"
// May need other stuff here

using namespace std;

namespace Hammer {

    TEST(PhaseSpaceTest, PS) {

        double com3 = compareVals(phaseSpaceN(1.2, {0.2, 0.3, 0.01}), phaseSpaceNBody(1.2, {0.2, 0.3, 0.01}));
        EXPECT_NEAR(com3, 1., 1e-3);

        double com4 = compareVals(phaseSpaceN(1.2, {0.2, 0.3, 0.01, 0.4}), phaseSpaceNBody(1.2, {0.2, 0.3, 0.01, 0.4}));
        EXPECT_NEAR(com4, 1., 1e-3);

        double com5 =
            compareVals(phaseSpaceN(1.2, {0.2, 0.3, 0.01, 0.2, 0.1}), phaseSpaceNBody(1.2, {0.2, 0.3, 0.01, 0.2, 0.1}));
        EXPECT_NEAR(com5, 1., 1e-3);
    }


} // namespace Hammer
