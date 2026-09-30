///
/// @file  TestAddAt.cc
/// @brief Tests for AddAt
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

    TEST(AddAtTest, SparseSparse) {
        auto t1 = makeEmptySparse({2, 2, 2, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        auto t2 =
            makeEmptySparse({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 0, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.2;
        t1->element({0, 0, 0, 0}) = 0.1i;
        t1->element({0, 1, 1, 1}) = 0.2;
        t1->element({0, 0, 1, 1}) = 0.2;
        t1->element({0, 1, 0, 1}) = 0.2;
        t1->element({1, 1, 1, 1}) = 0.5;
        auto res = addAt(std::move(t2), *t1, 2, 7);
        auto res2 = addAt(std::move(res), *t1, 2, 1);
        IContainer* tmp = res2.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 2, 1, 1}).real(), 0.);
    }

    TEST(AddAtTest, SparseVector) {
        auto t1 = makeEmptyVector({2, 2, 2, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        auto t2 =
            makeEmptySparse({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 0, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.2;
        t1->element({0, 0, 0, 0}) = 0.1i;
        t1->element({0, 1, 1, 1}) = 0.2;
        t1->element({0, 0, 1, 1}) = 0.2;
        t1->element({0, 1, 0, 1}) = 0.2;
        t1->element({1, 1, 1, 1}) = 0.5;
        auto res = addAt(std::move(t2), *t1, 2, 7);
        auto res2 = addAt(std::move(res), *t1, 2, 1);
        IContainer* tmp = res2.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 2, 1, 1}).real(), 0.);
    }

    TEST(AddAtTest, SparseVector2) {
        auto t1 = makeEmptyVector({2, 2, 2, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        auto t2 =
            makeEmptySparse({11, 2, 2, 2, 2}, {INTEGRATION_INDEX, SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 0, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.2;
        t1->element({0, 0, 0, 0}) = 0.1i;
        t1->element({0, 1, 1, 1}) = 0.2;
        t1->element({0, 0, 1, 1}) = 0.2;
        t1->element({0, 1, 0, 1}) = 0.2;
        t1->element({1, 1, 1, 1}) = 0.5;
        auto res = addAt(std::move(t2), *t1, 0, 7);
        auto res2 = addAt(std::move(res), *t1, 0, 1);
        IContainer* tmp = res2.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_DOUBLE_EQ(res2->element({7, 1, 0, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({7, 1, 0, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({7, 0, 0, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({7, 0, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({7, 0, 0, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({7, 0, 1, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({7, 1, 1, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 0, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 0, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 0, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 0, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({2, 1, 0, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({2, 1, 0, 1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({2, 0, 0, 0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({2, 0, 1, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({2, 0, 0, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({2, 0, 1, 0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({2, 1, 1, 1, 1}).real(), 0.);
    }

    TEST(AddAtTest, VectorSparse) {
        auto t1 = makeEmptySparse({2, 2, 2, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        auto t2 =
            makeEmptyVector({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 0, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.2;
        t1->element({0, 0, 0, 0}) = 0.1i;
        t1->element({0, 1, 1, 1}) = 0.2;
        t1->element({0, 0, 1, 1}) = 0.2;
        t1->element({0, 1, 0, 1}) = 0.2;
        t1->element({1, 1, 1, 1}) = 0.5;
        auto res = addAt(std::move(t2), *t1, 2, 7);
        auto res2 = addAt(std::move(res), *t1, 2, 1);
        IContainer* tmp = res2.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_NE(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 2, 1, 1}).real(), 0.);
    }

    TEST(AddAtTest, VectorVector) {
        auto t1 = makeEmptyVector({2, 2, 2, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        auto t2 =
            makeEmptyVector({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 0, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.2;
        t1->element({0, 0, 0, 0}) = 0.1i;
        t1->element({0, 1, 1, 1}) = 0.2;
        t1->element({0, 0, 1, 1}) = 0.2;
        t1->element({0, 1, 0, 1}) = 0.2;
        t1->element({1, 1, 1, 1}) = 0.5;
        auto res = addAt(std::move(t2), *t1, 2, 7);
        auto res2 = addAt(std::move(res), *t1, 2, 1);
        IContainer* tmp = res2.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_NE(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 0, 0}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 1, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 0, 1}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 1, 1}).real(), 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 2, 1, 1}).real(), 0.);
    }

    TEST(AddAtTest, SparseOuter) {
        auto t1 = makeEmptyVector({2, 2}, {SPIN_TAUP, SPIN_MUP});
        auto t2 =
            makeEmptySparse({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0}) = 0.1;
        t1->element({0, 0}) = 0.1i;
        t1->element({0, 1}) = 0.2;
        t1->element({1, 1}) = 0.5;
        auto t3 = makeOuterSquare(t1);
        auto res = addAt(std::move(t2), *t3, 2, 7);
        auto res2 = addAt(std::move(res), *t3, 2, 1);
        IContainer* tmp = res2.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 0, 0}).imag(), 0.1 * (-0.1));
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 0, 0}).real(), 0.1 * (-0.1) * -1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 0, 0}).imag(), 0.2 * (-0.1));
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 0, 0}).imag(), 0.5 * (-0.1));
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 1, 1}).real(), 0.1 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 1, 1}).imag(), 0.1 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 1, 1}).real(), 0.2 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 1, 1}).real(), 0.5 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 0, 1}).real(), 0.1 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 0, 1}).imag(), 0.1 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 0, 1}).real(), 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 0, 1}).real(), 0.5 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 1, 0}).real(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 1, 0}).imag(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 1, 0}).real(), 0.2 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 1, 0}).real(), 0.5 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 2, 1, 1}).real(), 0.);
    }


    TEST(AddAtTest, VectorOuter) {
        auto t1 = makeEmptyVector({2, 2}, {SPIN_TAUP, SPIN_MUP});
        auto t2 =
            makeEmptyVector({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0}) = 0.1;
        t1->element({0, 0}) = 0.1i;
        t1->element({0, 1}) = 0.2;
        t1->element({1, 1}) = 0.5;
        auto t3 = makeOuterSquare(t1);
        auto res = addAt(std::move(t2), *t3, 2, 7);
        auto res2 = addAt(std::move(res), *t3, 2, 1);
        IContainer* tmp = res2.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_NE(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 0, 0}).imag(), 0.1 * (-0.1));
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 0, 0}).real(), 0.1 * (-0.1) * -1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 0, 0}).imag(), 0.2 * (-0.1));
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 0, 0}).imag(), 0.5 * (-0.1));
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 7, 1, 1}).real(), 0.1 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 7, 1, 1}).imag(), 0.1 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 7, 1, 1}).real(), 0.2 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 7, 1, 1}).real(), 0.5 * 0.5);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 0, 1}).real(), 0.1 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 0, 1}).imag(), 0.1 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 0, 1}).real(), 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 0, 1}).real(), 0.5 * 0.2);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 1, 1, 0}).real(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 1, 1, 0}).imag(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 1, 1, 0}).real(), 0.2 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 1, 1, 0}).real(), 0.5 * 0.1);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 0, 2, 1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 0, 2, 1, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({0, 1, 2, 0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(res2->element({1, 1, 2, 1, 1}).real(), 0.);
    }


} // namespace Hammer::MultiDimensional
