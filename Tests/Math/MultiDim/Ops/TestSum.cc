///
/// @file  TestSum.cc
/// @brief Tests for Sum
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

    TEST(SumTest, SparseSparse) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        auto t2 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        auto res = sum(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.4);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 + 0.3);
    }

    TEST(SumTest, VectorVector) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        auto t2 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        auto res = sum(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.4);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 + 0.3);
    }

    TEST(SumTest, VectorSparse) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        auto t2 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        auto res = sum(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.4);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 + 0.3);
    }

    TEST(SumTest, SparseVector) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2;
        auto t2 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t2->element({1, 2, 0}) = 0.3;
        t2->element({1, 0, 0}) = 0.4;
        auto res = sum(std::move(t1), *t2);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 0}).real(), 0.4);
        EXPECT_DOUBLE_EQ(res->element({1, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 2, 0}).real(), 0.1 + 0.3);
    }

    TEST(SumTest, OuterOuter) {
        auto t = makeEmptySparse({2, 3}, {WILSON_BCENU, FF_BD});
        t->element({1, 2}) = 0.1;
        t->element({0, 1}) = 0.2;
        auto t2 = makeOuterSquare(t);
        auto t2b = t2->clone();
        auto t3 = makeEmptySparse({2, 3}, {WILSON_BCENU, FF_BD});
        t3->element({1, 2}) = 0.4;
        t3->element({0, 1}) = 0.7;
        t3->element({0, 0}) = 0.1;
        auto t4 = makeOuterSquare(t3);
        auto res = sum(std::move(t2), *t4);
        /// @todo fix this const crap
        const OuterContainer* o2 = static_cast<OuterContainer*>(t2b.get());
        const OuterContainer* o4 = static_cast<OuterContainer*>(t4.get());
        const OuterContainer* ores = static_cast<OuterContainer*>(res.get());
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        EXPECT_DOUBLE_EQ(ores->element({i1, i2, i3, i4}).real(),
                                         o2->element({i1, i2, i3, i4}).real() + o4->element({i1, i2, i3, i4}).real());
                        EXPECT_DOUBLE_EQ(ores->element({i1, i2, i3, i4}).imag(),
                                         o2->element({i1, i2, i3, i4}).imag() + o4->element({i1, i2, i3, i4}).imag());
                    }
                }
            }
        }
    }

    TEST(SumTest, OuterAny) {
        auto t = makeEmptySparse({2, 3}, {WILSON_BCENU, FF_BD});
        t->element({1, 2}) = 0.1;
        t->element({0, 1}) = 0.2;
        auto t2 = makeOuterSquare(t);
        auto t2b = t2->clone();
        auto t3 = makeEmptySparse({2, 3, 2, 3}, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC});
        t3->element({1, 2, 0, 0}) = 0.4 * 0.1;
        t3->element({1, 2, 0, 1}) = 0.4 * 0.7;
        t3->element({1, 2, 1, 2}) = 0.4 * 0.4;
        t3->element({0, 1, 0, 0}) = 0.7 * 0.1;
        t3->element({0, 1, 0, 1}) = 0.7 * 0.7;
        t3->element({0, 1, 1, 2}) = 0.7 * 0.4;
        t3->element({0, 0, 0, 0}) = 0.1 * 0.1;
        t3->element({0, 0, 0, 1}) = 0.1 * 0.7;
        t3->element({0, 0, 1, 2}) = 0.1 * 0.4;
        auto res = sum(std::move(t2), *t3);
        /// @todo fix this const crap
        const OuterContainer* o2 = static_cast<OuterContainer*>(t2b.get());
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        EXPECT_DOUBLE_EQ(res->element({i1, i2, i3, i4}).real(),
                                         o2->element({i1, i2, i3, i4}).real() + t3->element({i1, i2, i3, i4}).real());
                        EXPECT_DOUBLE_EQ(res->element({i1, i2, i3, i4}).imag(),
                                         o2->element({i1, i2, i3, i4}).imag() + t3->element({i1, i2, i3, i4}).imag());
                    }
                }
            }
        }
    }

} // namespace Hammer::MultiDimensional
