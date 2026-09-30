///
/// @file  TestMathUtils.cc
/// @brief Tests for MathUtils
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include "Hammer/Math/Utils.hh"

#include <complex>

using namespace std;

namespace Hammer {

    TEST(MathUtilsTest, compareValsDoubleNonZero) {
        EXPECT_DOUBLE_EQ(compareVals(2.0, 4.0), 0.5);
        EXPECT_DOUBLE_EQ(compareVals(3.0, 1.0), 3.0);
    }

    TEST(MathUtilsTest, compareValsDoubleZero) {
        EXPECT_DOUBLE_EQ(compareVals(0.0, 0.0), 1.0);
        EXPECT_DOUBLE_EQ(compareVals(3.0, 0.0), 4.0);
    }

    TEST(MathUtilsTest, compareValsComplexNonZero) {
        complex<double> num{2., 1.};
        complex<double> den{4., 2.};
        complex<double> result = compareVals(num, den);
        EXPECT_DOUBLE_EQ(result.real(), 0.5);
        EXPECT_DOUBLE_EQ(result.imag(), 0.0);
    }

    TEST(MathUtilsTest, compareValsComplexZero) {
        complex<double> zero{0., 0.};
        complex<double> result = compareVals(zero, zero);
        EXPECT_DOUBLE_EQ(result.real(), 1.0);
        EXPECT_DOUBLE_EQ(result.imag(), 0.0);
    }

} // namespace Hammer
