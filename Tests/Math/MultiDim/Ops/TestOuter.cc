///
/// @file  TestOuter.cc
/// @brief Tests for Outer
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

    TEST(OuterTest, Sparse) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2 + 0.3i;
        auto res = calcSquare(std::move(t1));
        EXPECT_EQ(res->rank(), 6);
        EXPECT_EQ(res->dims()[3], 2);
        EXPECT_EQ(res->dims()[4], 3);
        EXPECT_EQ(res->dims()[5], 1);
        EXPECT_EQ(res->labels()[3], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[4], FF_BD_HC);
        EXPECT_EQ(res->labels()[5], SPIN_DSTAR_HC);
        const OuterContainer* ores = static_cast<OuterContainer*>(res.get());
        EXPECT_DOUBLE_EQ(ores->element({0, 0, 0, 0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(ores->element({0, 1, 0, 0, 1, 0}).real(), 0.2 * 0.2 + 0.3 * 0.3);
        EXPECT_DOUBLE_EQ(ores->element({0, 2, 0, 0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(ores->element({1, 2, 0, 1, 2, 0}).real(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(ores->element({1, 2, 0, 0, 1, 0}).real(), 0.1 * 0.2);
        EXPECT_DOUBLE_EQ(ores->element({1, 2, 0, 0, 1, 0}).imag(), -0.1 * 0.3);
    }

    TEST(OuterTest, Vector) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2 + 0.3i;
        auto res = calcSquare(std::move(t1));
        const OuterContainer* ores = static_cast<OuterContainer*>(res.get());
        EXPECT_DOUBLE_EQ(ores->element({0, 0, 0, 0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(ores->element({0, 1, 0, 0, 1, 0}).real(), 0.2 * 0.2 + 0.3 * 0.3);
        EXPECT_DOUBLE_EQ(ores->element({0, 2, 0, 0, 2, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(ores->element({1, 2, 0, 1, 2, 0}).real(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(ores->element({1, 2, 0, 0, 1, 0}).real(), 0.1 * 0.2);
        EXPECT_DOUBLE_EQ(ores->element({1, 2, 0, 0, 1, 0}).imag(), -0.1 * 0.3);
    }

    TEST(OuterTest, Outer) {
        auto t = makeEmptySparse({2}, {WILSON_BCENU});
        t->element({1}) = 0.1;
        t->element({0}) = 0.2i;
        auto t2 = makeOuterSquare(t);
        auto res = calcSquare(std::move(t2));
        EXPECT_EQ(res->labels()[0], WILSON_BCENU);
        EXPECT_EQ(res->labels()[1], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[2], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[3], WILSON_BCENU);
        /// @todo fix this const crap
        const OuterContainer* ores = static_cast<OuterContainer*>(res.get());
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 2; ++i4) {
                        EXPECT_DOUBLE_EQ(
                            ores->element({i1, i2, i3, i4}).real(),
                            (t->element({i1}) * conj(t->element({i2}) * t->element({i3})) * t->element({i4})).real());
                        EXPECT_DOUBLE_EQ(
                            ores->element({i1, i2, i3, i4}).imag(),
                            (t->element({i1}) * conj(t->element({i2}) * t->element({i3})) * t->element({i4})).imag());
                    }
                }
            }
        }
    }

    TEST(OuterTest, Outer2) {
        auto t = makeEmptySparse({2}, {WILSON_BCENU});
        t->element({1}) = 0.1;
        t->element({0}) = 0.2i;
        auto t2 = makeEmptySparse({2}, {WILSON_BCENU});
        t2->element({1}) = 0.7;
        t2->element({0}) = 0.9i;
        auto t3 = makeOuterSquare(t);
        auto t4 = makeOuterSquare(t2);
        auto tmp = sum(std::move(t3), *t4);
        auto tmpb = tmp->clone();
        auto res = calcSquare(std::move(tmp));
        EXPECT_EQ(res->labels()[0], WILSON_BCENU);
        EXPECT_EQ(res->labels()[1], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[2], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[3], WILSON_BCENU);
        /// @todo fix this const crap
        const OuterContainer* ores = static_cast<OuterContainer*>(res.get());
        const OuterContainer* otmp = static_cast<OuterContainer*>(tmpb.get());
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 2; ++i4) {
                        EXPECT_DOUBLE_EQ(ores->element({i1, i2, i3, i4}).real(),
                                         (otmp->element({i1, i2}) * conj(otmp->element({i3, i4}))).real());
                        EXPECT_DOUBLE_EQ(ores->element({i1, i2, i3, i4}).imag(),
                                         (otmp->element({i1, i2}) * conj(otmp->element({i3, i4}))).imag());
                    }
                }
            }
        }
    }

} // namespace Hammer::MultiDimensional
