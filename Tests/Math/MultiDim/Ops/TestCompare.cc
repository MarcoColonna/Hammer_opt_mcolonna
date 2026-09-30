///
/// @file  TestCompare.cc
/// @brief Tests for Compare
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/Operations.hh"
#include "Hammer/Exceptions.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    TEST(CompareTest, SparseVector) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2 + 0.3i;
        auto t2 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t2->element({1, 2, 0}) = 0.1;
        t2->element({0, 1, 0}) = 0.2 + 0.3i;
        auto t3 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t3->element({1, 2, 0}) = 0.1;
        t3->element({0, 1, 0}) = 0.2 + 0.4i;
        EXPECT_TRUE(*t1 == *t2);
        EXPECT_FALSE(*t1 == *t3);
        EXPECT_TRUE(*t2 == *t1);
        EXPECT_FALSE(*t3 == *t1);
    }

    TEST(CompareTest, SparseSparse) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2 + 0.3i;
        auto t2 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t2->element({1, 2, 0}) = 0.1;
        t2->element({0, 1, 0}) = 0.2 + 0.3i;
        auto t3 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t3->element({1, 2, 0}) = 0.1;
        t3->element({0, 1, 0}) = 0.2 + 0.4i;
        EXPECT_TRUE(*t1 == *t2);
        EXPECT_FALSE(*t1 == *t3);
    }

    TEST(CompareTest, VectorVector) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2 + 0.3i;
        auto t2 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t2->element({1, 2, 0}) = 0.1;
        t2->element({0, 1, 0}) = 0.2 + 0.3i;
        auto t3 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t3->element({1, 2, 0}) = 0.1;
        t3->element({0, 1, 0}) = 0.2 + 0.4i;
        EXPECT_TRUE(*t1 == *t2);
        EXPECT_FALSE(*t1 == *t3);
    }

    TEST(CompareTest, VectorOuter) {
        auto t1 = makeEmptyVector({2, 2}, {WILSON_BCENU, WILSON_BCENU_HC});
        t1->element({0, 0}) = (0.2i) * (-0.2i) + (0.9i) * (-0.9i);
        t1->element({0, 1}) = (0.2i) * 0.1 + (0.9i) * 0.7;
        t1->element({1, 0}) = 0.1 * (-0.2i) + 0.7 * (-0.9i);
        t1->element({1, 1}) = 0.1 * 0.1 + 0.7 * 0.7;
        auto t = makeEmptyVector({2}, {WILSON_BCENU});
        t->element({1}) = 0.1;
        t->element({0}) = 0.2i;
        auto t2 = makeEmptyVector({2}, {WILSON_BCENU});
        t2->element({1}) = 0.7;
        t2->element({0}) = 0.9i;
        auto t3 = makeOuterSquare(t);
        auto t4 = makeOuterSquare(t2);
        auto tmp = sum(std::move(t3), *t4);
        EXPECT_TRUE(*t1 == *tmp);
    }

} // namespace Hammer::MultiDimensional
