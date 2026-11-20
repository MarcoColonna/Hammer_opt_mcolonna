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
        auto b = BruteForceIteratorRange(t.dims());
        auto b2 = BruteForceIteratorRange(t.dims(), {2, 1, 1});
        auto b3 = BruteForceIteratorRange(t.dims(), {1, 1, 1});
        auto b4 = BruteForceIteratorRange(t.dims(), {2, 2, 0});
        IndexList s{0, 0, 0};
        IndexList s2{0, 1, 0};
        IndexList s3{1, 1, 0};
        IndexList s4{0, 2, 0};

        EXPECT_EQ(*b.begin(), s);
        EXPECT_EQ(*b2.begin(), s2);
        EXPECT_EQ(*b3.begin(), s3);
        EXPECT_EQ(*b4.begin(), s4);
    }

    TEST(BruteForceTest, FullyFixedCoordinates) {
        const IndexList dimensions{3, 2};
        for (const auto& fixed : {IndexList{0, 0}, IndexList{2, 1}}) {
            const BruteForceIteratorRange range{dimensions, fixed};
            auto it = range.begin();
            ASSERT_NE(it, range.end());
            EXPECT_EQ(*it, fixed);
            EXPECT_EQ(++it, range.end());
        }
    }

    TEST(BruteForceTest, Increment) {
        SequentialIndexing t{{2, 3, 1}};
        auto b = BruteForceIteratorRange(t.dims());
        auto b2 = BruteForceIteratorRange(t.dims(), {2, 1, 1});
        auto b3 = BruteForceIteratorRange(t.dims(), {1, 1, 1});
        auto b4 = BruteForceIteratorRange(t.dims(), {2, 2, 0});
        auto b5 = BruteForceIteratorRange(t.dims(), {1, 3, 0});
        EXPECT_TRUE(b.begin() < b.end());
        EXPECT_EQ(distance(b.begin(), b.end()), 6);
        EXPECT_EQ(distance(b2.begin(), b2.end()), 2);
        EXPECT_EQ(distance(b3.begin(), b3.end()), 1);
        EXPECT_EQ(distance(b4.begin(), b4.end()), 2);
        EXPECT_EQ(distance(b5.begin(), b5.end()), 3);
    }


} // namespace Hammer::MultiDimensional
