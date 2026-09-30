///
/// @file  TestBruteForceIterator.cc
/// @brief Tests for BruteForceIterator
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/SequentialIndexing.hh"
#include "Hammer/Math/MultiDim/BruteForceIterator.hh"
#include "Hammer/Exceptions.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    TEST(BruteForceTest, ConstructionAccess) {
        SequentialIndexing t{{2, 3, 1}};
        auto b = BruteForceIterator(t.dims());
        auto b2 = BruteForceIterator(t.dims(), {2, 1, 1});
        auto b3 = BruteForceIterator(t.dims(), {1, 1, 1});
        auto b4 = BruteForceIterator(t.dims(), {2, 2, 0});
        IndexList s{0, 0, 0};
        IndexList s2{0, 1, 0};
        IndexList s3{1, 1, 0};
        IndexList s4{0, 2, 0};

        EXPECT_EQ(*b, s);
        EXPECT_EQ(*b2, s2);
        EXPECT_EQ(*b3, s3);
        EXPECT_EQ(*b4, s4);
    }

    TEST(BruteForceTest, Increment) {
        SequentialIndexing t{{2, 3, 1}};
        auto b = BruteForceIterator(t.dims());
        auto b2 = BruteForceIterator(t.dims(), {2, 1, 1});
        auto b3 = BruteForceIterator(t.dims(), {1, 1, 1});
        auto b4 = BruteForceIterator(t.dims(), {2, 2, 0});
        auto b5 = BruteForceIterator(t.dims(), {1, 3, 0});
        EXPECT_THROW(std::ignore = distance(b.end(), b.begin()), RangeError);
        EXPECT_EQ(distance(b.begin(), b.end()), 6);
        EXPECT_EQ(distance(b2.begin(), b2.end()), 2);
        EXPECT_EQ(distance(b3.begin(), b3.end()), 1);
        EXPECT_EQ(distance(b4.begin(), b4.end()), 2);
        EXPECT_EQ(distance(b5.begin(), b5.end()), 3);
    }


} // namespace Hammer::MultiDimensional
