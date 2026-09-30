///
/// @file  TestScalarContainer.cc
/// @brief Tests for ScalarContainer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Exceptions.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    TEST(ScalarTest, ConstructionAccess) {
        auto t = makeEmptyScalar();
        EXPECT_EQ(t->rank(), 0);
    }

    TEST(ScalarTest, ElementAccess) {
        auto t = makeEmptyScalar();
        t->element({}) = 0.1;
        EXPECT_THROW(t->element({1, 2, 1}), RangeError);
        EXPECT_DOUBLE_EQ(t->element({}).real(), 0.1);
    }

    TEST(ScalarTest, Comparison) {
        auto t1 = makeEmptyScalar();
        auto t2 = makeEmptyScalar();
        EXPECT_TRUE(t1->compare(*t2));
        t1->element({}) = 0.3;
        EXPECT_FALSE(t1->compare(*t2));
    }

    TEST(ScalarTest, ComponentMultiply) {
        auto t = makeEmptyScalar();
        t->element({}) = 0.1;
        EXPECT_DOUBLE_EQ(t->element({}).real(), 0.1);
        (*t) *= 3.;
        EXPECT_DOUBLE_EQ(t->element({}).real(), 0.3);
    }

    TEST(ScalarTest, Cloning) {
        auto t = makeEmptyScalar();
        t->element({}) = 0.1;
        auto t2 = t->clone();
        EXPECT_EQ(t->rank(), t2->rank());
        EXPECT_DOUBLE_EQ(t->element({}).real(), t2->element({}).real());
        EXPECT_DOUBLE_EQ(t->element({}).imag(), t2->element({}).imag());
    }

    TEST(ScalarTest, Conjugate) {
        auto t = makeEmptyScalar();
        t->element({}) = 1.0 + 0.2i;
        t->conjugate();
        EXPECT_DOUBLE_EQ(t->element({}).real(), 1.0);
        EXPECT_DOUBLE_EQ(t->element({}).imag(), -0.2);
    }

    TEST(ScalarTest, QueryMethods) {
        auto t = makeEmptyScalar();
        EXPECT_EQ(t->numValues(), 1ul);
        EXPECT_EQ(t->dataSize(), sizeof(complex<double>));
        EXPECT_EQ(t->entrySize(), sizeof(complex<double>));
    }

    TEST(ScalarTest, LabelToIndex) {
        auto t = makeEmptyScalar();
        EXPECT_EQ(t->labelToIndex(WILSON_BCENU), 0ul);
    }

    TEST(ScalarTest, LabelPairQueries) {
        auto t1 = makeEmptyScalar();
        auto t2 = makeEmptyScalar();
        UniqueLabelsList ulist{};
        auto pairs = t1->getSameLabelPairs(*t2, ulist);
        EXPECT_TRUE(pairs.empty());
        auto spinPairs = t1->getSpinLabelPairs();
        EXPECT_TRUE(spinPairs.empty());
    }

    TEST(ScalarTest, IsSameShape) {
        auto t1 = makeEmptyScalar();
        auto t2 = makeEmptyScalar();
        EXPECT_TRUE(t1->isSameShape(*t2));
        auto t3 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        EXPECT_FALSE(t1->isSameShape(*t3));
    }

    TEST(ScalarTest, CanAddAt) {
        auto t = makeEmptyScalar();
        auto tsub = makeEmptyScalar();
        EXPECT_FALSE(t->canAddAt(*tsub, WILSON_BCENU, 0));
    }

    TEST(ScalarTest, ElementRangeIterator) {
        auto t = makeEmptyScalar();
        t->element({}) = 0.5 + 0.3i;
        IndexList empty{};
        EXPECT_DOUBLE_EQ(t->element(empty.begin(), empty.end()).real(), 0.5);
        EXPECT_DOUBLE_EQ(t->element(empty.begin(), empty.end()).imag(), 0.3);
        IndexList nonempty{1};
        EXPECT_THROW(t->element(nonempty.begin(), nonempty.end()), RangeError);
    }

    TEST(ScalarTest, ElementRangeIteratorConst) {
        auto t = makeEmptyScalar();
        t->element({}) = 0.7 - 0.1i;
        const IContainer* ct = t.get();
        IndexList empty{};
        EXPECT_DOUBLE_EQ(ct->element(empty.begin(), empty.end()).real(), 0.7);
        EXPECT_DOUBLE_EQ(ct->element(empty.begin(), empty.end()).imag(), -0.1);
        IndexList nonempty{0};
        EXPECT_THROW({ std::ignore = ct->element(nonempty.begin(), nonempty.end()); }, RangeError);
    }

    TEST(ScalarTest, ComparisonWithNonScalar) {
        auto t1 = makeEmptyScalar();
        auto t3 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        EXPECT_FALSE(t1->compare(*t3));
    }

} // namespace Hammer::MultiDimensional
