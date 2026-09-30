///
/// @file  TestDifferentiator.cc
/// @brief Tests for Differentiator
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/Differentiator.hh"
#include "gtest/gtest.h"
// May need other stuff here

using namespace std;

namespace Hammer {

    TEST(DifferentiatorTest, Center) {

        std::function<double(double)> f = [](double x) -> double { return 4. * x + 5. * x * x - 2. + 8. * x * x * x; };
        EXPECT_NEAR(der<1>(f, 0.), 4., 1.e-6);
        EXPECT_NEAR(der<1>(f, 1.), 38., 1.e-6);
        EXPECT_NEAR(der<1>(f, -1.), 18., 1.e-6);
        EXPECT_NEAR(der<2>(f, 0.), 10., 1.e-3);
        EXPECT_NEAR(der<2>(f, 1.), 58., 1.e-3);
        EXPECT_NEAR(der<2>(f, -1.), -38., 1.e-3);
        EXPECT_THROW(der<3>(f, 1.), Error);
    }

    TEST(DifferentiatorTest, Right) {

        std::function<double(double)> f = [](double x) -> double { return 4. * x + 5. * x * x - 2. + 8. * x * x * x; };
        EXPECT_NEAR(der<1>(f, 0., {0., 3.}), 4., 1.e-6);
        EXPECT_NEAR(der<1>(f, 1., {1., 3.}), 38., 1.e-6);
        EXPECT_NEAR(der<1>(f, -1., {-1., 3.}), 18., 1.e-6);
        EXPECT_NEAR(der<2>(f, 0., {0., 3.}), 10., 1.e-3);
        EXPECT_NEAR(der<2>(f, 1., {1., 3.}), 58., 1.e-3);
        EXPECT_NEAR(der<2>(f, -1., {-1., 3.}), -38., 1.e-3);
        EXPECT_THROW(der<3>(f, 1., {1., 3.}), Error);
    }

    TEST(DifferentiatorTest, Left) {

        std::function<double(double)> f = [](double x) -> double { return 4. * x + 5. * x * x - 2. + 8. * x * x * x; };
        EXPECT_NEAR(der<1>(f, 0., {-3., 0.}), 4., 1.e-6);
        EXPECT_NEAR(der<1>(f, 1., {-3., 1.}), 38., 1.e-6);
        EXPECT_NEAR(der<1>(f, -1., {-3., -1.}), 18., 1.e-6);
        EXPECT_NEAR(der<2>(f, 0., {-3., 0.}), 10., 1.e-3);
        EXPECT_NEAR(der<2>(f, 1., {-3., 1.}), 58., 1.e-3);
        EXPECT_NEAR(der<2>(f, -1., {-3., -1.}), -38., 1.e-3);
        EXPECT_THROW(der<3>(f, 1., {1., 3.}), Error);
    }

} // namespace Hammer
