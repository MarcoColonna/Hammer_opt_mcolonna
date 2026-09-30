///
/// @file  TestUnits.cc
/// @brief Tests for Units
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include "Hammer/Math/Units.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer {

    TEST(UnitsTest, Conversions) {
        auto& u = Units::instance();
        EXPECT_DOUBLE_EQ(u.getUnitsRescalingToMC("GeV", "GeV"), 1.0);
        EXPECT_DOUBLE_EQ(u.getUnitsRescalingToMC("GeV", "MeV"), 1e-3);
        EXPECT_DOUBLE_EQ(u.getUnitsRescalingToMC("MeV", "GeV"), 1e3);
        EXPECT_DOUBLE_EQ(u.getUnitsRescalingToMC("MeV", "MeV"), 1.0);
        EXPECT_DOUBLE_EQ(u.getUnitsRescalingToMC("gev", "mev"), 1e-3);
    }

    TEST(UnitsTest, BadMcUnits) {
        auto& u = Units::instance();
        EXPECT_THROW(std::ignore = u.getUnitsRescalingToMC("Parsec", "GeV"), Error);
    }

    TEST(UnitsTest, BadLocalUnits) {
        auto& u = Units::instance();
        EXPECT_THROW(std::ignore = u.getUnitsRescalingToMC("GeV", "Parsec"), Error);
    }

} // namespace Hammer
