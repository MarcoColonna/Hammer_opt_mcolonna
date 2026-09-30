///
/// @file  TestVectorContainer.cc
/// @brief Tests for VectorContainer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    TEST(VectorTest, ConstructionAccess) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        EXPECT_EQ(t->rank(), 3);
        EXPECT_EQ(t->dims()[0], 2);
        EXPECT_EQ(t->dims()[1], 3);
        EXPECT_EQ(t->dims()[2], 1);
    }

    TEST(VectorTest, ElementAccess) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2;
#ifndef NDEBUG
        ASSERT_DEATH(t->element({1, 2, 1}), "");
#endif
        EXPECT_DOUBLE_EQ(t->element({1, 2, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(t->element({0, 1, 0}).real(), 0.2);
    }

    TEST(VectorTest, Iterators) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        unique_ptr<VectorContainer> t2{static_cast<VectorContainer*>(t.release())};
        auto it1 = t2->begin();
        *it1 = 0.1;
        auto it2 = t2->end();
        --it2;
        *it2 = 0.2;
        EXPECT_DOUBLE_EQ(t2->element({0, 0, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(t2->element({1, 2, 0}).real(), 0.2);
    }

    TEST(VectorTest, BracketAccess) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        unique_ptr<VectorContainer> t2{static_cast<VectorContainer*>(t.release())};
        (*t2)[1 * 3 + 2] = 0.1;
        (*t2)[0 * 3 + 1] = 0.2;
        EXPECT_DOUBLE_EQ(t2->element({1, 2, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(t2->element({0, 1, 0}).real(), 0.2);
    }

    TEST(VectorTest, ValueAccess) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        unique_ptr<VectorContainer> t2{static_cast<VectorContainer*>(t.release())};
        t2->setValue({1, 2, 0}, 0.1);
        t2->setValue({0, 1, 0}, 0.2);
        EXPECT_DOUBLE_EQ(t2->value({1, 2, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(t2->value({0, 1, 0}).real(), 0.2);
    }

    TEST(VectorTest, Comparison) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        auto t2 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        auto t3 = makeEmptyVector({4, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        auto t4 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD_HC, INTEGRATION_INDEX});
        auto t5 = makeEmptyScalar();
        EXPECT_TRUE(t1->compare(*t2));
        EXPECT_FALSE(t1->compare(*t3));
        EXPECT_FALSE(t1->compare(*t4));
        EXPECT_FALSE(t1->compare(*t5));
        t1->element({0, 1, 0}) = 0.3;
        EXPECT_FALSE(t1->compare(*t2));
    }

    TEST(VectorTest, ComponentMultiply) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2;
        EXPECT_DOUBLE_EQ(t->element({1, 2, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(t->element({0, 1, 0}).real(), 0.2);
        EXPECT_DOUBLE_EQ(t->element({0, 0, 0}).real(), 0.0);
        (*t) *= 3.;
        EXPECT_DOUBLE_EQ(t->element({1, 2, 0}).real(), 0.3);
        EXPECT_DOUBLE_EQ(t->element({0, 1, 0}).real(), 0.6);
        EXPECT_DOUBLE_EQ(t->element({0, 0, 0}).real(), 0.0);
    }

    TEST(VectorTest, Cloning) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2;
        auto t2 = t->clone();
        EXPECT_EQ(t->rank(), t2->rank());
        EXPECT_EQ(t->dims()[0], t2->dims()[0]);
        EXPECT_EQ(t->dims()[1], t2->dims()[1]);
        EXPECT_EQ(t->dims()[2], t2->dims()[2]);
        EXPECT_EQ(t->labels()[0], t2->labels()[0]);
        EXPECT_EQ(t->labels()[1], t2->labels()[1]);
        EXPECT_EQ(t->labels()[2], t2->labels()[2]);
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                EXPECT_DOUBLE_EQ(t->element({i1, i2, 0}).real(), t2->element({i1, i2, 0}).real());
                EXPECT_DOUBLE_EQ(t->element({i1, i2, 0}).imag(), t2->element({i1, i2, 0}).imag());
            }
        }
    }

    TEST(VectorTest, Conjugate) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_TAUM});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2i;
        t->conjugate();
        EXPECT_DOUBLE_EQ(t->element({1, 2, 0}).real(), 0.1);
        EXPECT_DOUBLE_EQ(t->element({0, 1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(t->element({0, 1, 0}).imag(), -0.2);
        EXPECT_EQ(t->labels()[0], WILSON_BCENU_HC);
        EXPECT_EQ(t->labels()[1], FF_BD_HC);
        EXPECT_EQ(t->labels()[2], SPIN_TAUM_HC);
    }

    TEST(VectorTest, ItSequential) {
        auto t = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_TAUM});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2i;
        auto* sp = static_cast<VectorContainer*>(t.get());
        auto itb = sp->firstNonZero();
        auto ite = sp->endNonZero();
        EXPECT_EQ(itb->position(), 1);
        EXPECT_EQ(itb->value().imag(), 0.2);
        EXPECT_EQ(itb->value().real(), 0.);
        itb->next();
        EXPECT_EQ(itb->position(), 1 * 3 + 2);
        EXPECT_EQ(itb->value().imag(), 0.);
        EXPECT_EQ(itb->value().real(), 0.1);
        itb->next();
        EXPECT_EQ(*itb, *ite);
        EXPECT_FALSE(itb->isAligned());
    }


} // namespace Hammer::MultiDimensional
