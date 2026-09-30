///
/// @file  TestMultiply.cc
/// @brief Tests for Multiply
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

    TEST(MultiplyTest, SparseSparse) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        t1->element({0, 0, 0}) = 0.7;
        auto t2 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        t2->element({0, 1, 0}) = 0.4;
        auto res = elementMultiply(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 * 0.3);
    }

    TEST(MultiplyTest, VectorVector) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        t1->element({0, 0, 0}) = 0.7;
        auto t2 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        t2->element({0, 1, 0}) = 0.4;
        auto res = elementMultiply(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 * 0.3);
    }

    TEST(MultiplyTest, VectorSparse) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        t1->element({0, 0, 0}) = 0.7;
        auto t2 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        t2->element({0, 1, 0}) = 0.4;
        auto res = elementMultiply(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 * 0.3);
    }

    TEST(MultiplyTest, SparseVector) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        t1->element({0, 0, 0}) = 0.7;
        auto t2 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        t2->element({0, 1, 0}) = 0.4;
        auto res = elementMultiply(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 * 0.3);
    }

    TEST(MultiplyTest, OuterAny) {
        auto t = makeEmptySparse({2, 3}, {WILSON_BCENU, FF_BD});
        t->element({1, 2}) = 0.1;
        t->element({0, 1}) = 0.2;
        auto t2 = makeOuterSquare(t);
        auto t2b = makeOuterSquare(t);
        auto t2c = makeOuterSquare(t);
        auto t3 = makeEmptySparse({2, 3, 2, 3}, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC});
        t3->element({1, 2, 1, 2}) = 0.1 * 0.1;
        t3->element({1, 2, 0, 1}) = 0.2 * 0.1;
        t3->element({0, 1, 1, 2}) = 0.2 * 0.1;
        t3->element({0, 1, 0, 1}) = 0.2 * 0.2;
        auto t4 = makeEmptySparse({2}, {WILSON_BCENU});
        auto t5 = makeEmptySparse({3, 2, 3}, {FF_BD, WILSON_BCENU_HC, FF_BD_HC});
        vector<TensorData> vec;
        vec.push_back(std::move(t4));
        vec.push_back(std::move(t5));
        auto t6 = combineTensors(vec);
        auto t7 = makeEmptySparse({2, 3}, {WILSON_BCENU, FF_BD});
        t7->element({1, 2}) = 0.4;
        t7->element({0, 1}) = 0.7;
        t7->element({0, 0}) = 0.1;
        auto t8 = makeOuterSquare(t7);
        auto res = elementMultiply(std::move(t2), *t3);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 1, 2}).real(), pow(t3->element({1, 2, 1, 2}).real(), 2.0));
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0, 1}).real(), pow(t3->element({1, 2, 0, 1}).real(), 2.0));
        EXPECT_DOUBLE_EQ(res->element({0, 1, 1, 2}).real(), pow(t3->element({0, 1, 1, 2}).real(), 2.0));
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0, 1}).real(), pow(t3->element({0, 1, 0, 1}).real(), 2.0));
        EXPECT_THROW(elementMultiply(std::move(t2b), *t6), Error);
        EXPECT_THROW(elementMultiply(std::move(t2c), *t8), Error);
    }

} // namespace Hammer::MultiDimensional
