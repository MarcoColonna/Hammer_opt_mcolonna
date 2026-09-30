///
/// @file  TestTrace.cc
/// @brief Tests for Trace
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

    TEST(TraceTest, SparseAll) {
        auto t1 = makeEmptySparse({2, 2, 2, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 0, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.2;
        t1->element({0, 0, 0, 0}) = 0.1i;
        t1->element({0, 1, 1, 1}) = 0.2;
        t1->element({0, 0, 1, 1}) = 0.2;
        t1->element({0, 1, 0, 1}) = 0.2;
        t1->element({1, 1, 1, 1}) = 0.5;
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 0);
        EXPECT_EQ(pos[1].first, 1);
        EXPECT_EQ(pos[0].second, 2);
        EXPECT_EQ(pos[1].second, 3);
        auto res = calcTrace(std::move(t1), pos);
        IContainer* tmp = res.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element().real(), 0.2 + 0.2 + 0.5);
        EXPECT_DOUBLE_EQ(res->element().imag(), 0.1);
    }

    TEST(TraceTest, Sparse) {
        auto t1 =
            makeEmptySparse({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 4, 0, 0}) = 0.1;
        t1->element({1, 0, 4, 1, 0}) = 0.2;
        t1->element({0, 0, 4, 0, 0}) = 0.1i;
        t1->element({0, 1, 4, 1, 1}) = 0.2;
        t1->element({0, 0, 4, 1, 1}) = 0.2;
        t1->element({0, 1, 4, 0, 1}) = 0.2;
        t1->element({1, 1, 4, 1, 1}) = 0.5;
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 0);
        EXPECT_EQ(pos[1].first, 1);
        EXPECT_EQ(pos[0].second, 3);
        EXPECT_EQ(pos[1].second, 4);
        auto res = calcTrace(std::move(t1), pos);
        IContainer* tmp = res.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_EQ(res->rank(), 1);
        EXPECT_DOUBLE_EQ(res->element({4}).real(), 0.2 + 0.2 + 0.5);
        EXPECT_DOUBLE_EQ(res->element({4}).imag(), 0.1);
    }

    TEST(TraceTest, Sparse2) {
        auto t1 = makeEmptySparse({2, 2, 11, 2, 2, 11},
                                  {SPIN_TAUP, SPIN_MUP, WILSON_BCENU, SPIN_TAUP_HC, SPIN_MUP_HC, WILSON_BCENU_HC});
        t1->element({1, 0, 4, 0, 0, 4}) = 0.1;
        t1->element({1, 0, 4, 1, 0, 4}) = 0.2;
        t1->element({0, 0, 4, 0, 0, 4}) = 0.1i;
        t1->element({0, 1, 4, 1, 1, 4}) = 0.2;
        t1->element({0, 0, 4, 1, 1, 4}) = 0.2;
        t1->element({0, 1, 3, 0, 1, 2}) = 0.2;
        t1->element({1, 1, 2, 1, 1, 2}) = 0.5;
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 0);
        EXPECT_EQ(pos[1].first, 1);
        EXPECT_EQ(pos[0].second, 3);
        EXPECT_EQ(pos[1].second, 4);
        auto res = calcTrace(std::move(t1), pos);
        IContainer* tmp = res.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_EQ(res->rank(), 2);
        EXPECT_DOUBLE_EQ(res->element({4, 4}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({4, 4}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res->element({3, 2}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).real(), 0.5);
    }

    TEST(TraceTest, VectorAll) {
        auto t1 = makeEmptyVector({2, 2, 2, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 0, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.2;
        t1->element({0, 0, 0, 0}) = 0.1i;
        t1->element({0, 1, 1, 1}) = 0.2;
        t1->element({0, 0, 1, 1}) = 0.2;
        t1->element({0, 1, 0, 1}) = 0.2;
        t1->element({1, 1, 1, 1}) = 0.5;
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 0);
        EXPECT_EQ(pos[1].first, 1);
        EXPECT_EQ(pos[0].second, 2);
        EXPECT_EQ(pos[1].second, 3);
        auto res = calcTrace(std::move(t1), pos);
        IContainer* tmp = res.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_EQ(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_NE(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element().real(), 0.2 + 0.2 + 0.5);
        EXPECT_DOUBLE_EQ(res->element().imag(), 0.1);
    }

    TEST(TraceTest, Vector) {
        auto t1 =
            makeEmptyVector({2, 2, 11, 2, 2}, {SPIN_TAUP, SPIN_MUP, INTEGRATION_INDEX, SPIN_TAUP_HC, SPIN_MUP_HC});
        t1->element({1, 0, 4, 0, 0}) = 0.1;
        t1->element({1, 0, 4, 1, 0}) = 0.2;
        t1->element({0, 0, 4, 0, 0}) = 0.1i;
        t1->element({0, 1, 4, 1, 1}) = 0.2;
        t1->element({0, 0, 4, 1, 1}) = 0.2;
        t1->element({0, 1, 4, 0, 1}) = 0.2;
        t1->element({1, 1, 4, 1, 1}) = 0.5;
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 0);
        EXPECT_EQ(pos[1].first, 1);
        EXPECT_EQ(pos[0].second, 3);
        EXPECT_EQ(pos[1].second, 4);
        auto res = calcTrace(std::move(t1), pos);
        IContainer* tmp = res.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_NE(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_EQ(res->rank(), 1);
        EXPECT_DOUBLE_EQ(res->element({4}).real(), 0.2 + 0.2 + 0.5);
        EXPECT_DOUBLE_EQ(res->element({4}).imag(), 0.1);
    }

    TEST(TraceTest, Vector2) {
        auto t1 = makeEmptyVector({2, 2, 11, 2, 2, 11},
                                  {SPIN_TAUP, SPIN_MUP, WILSON_BCENU, SPIN_TAUP_HC, SPIN_MUP_HC, WILSON_BCENU_HC});
        t1->element({1, 0, 4, 0, 0, 4}) = 0.1;
        t1->element({1, 0, 4, 1, 0, 4}) = 0.2;
        t1->element({0, 0, 4, 0, 0, 4}) = 0.1i;
        t1->element({0, 1, 4, 1, 1, 4}) = 0.2;
        t1->element({0, 0, 4, 1, 1, 4}) = 0.2;
        t1->element({0, 1, 3, 0, 1, 2}) = 0.2;
        t1->element({1, 1, 2, 1, 1, 2}) = 0.5;
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 0);
        EXPECT_EQ(pos[1].first, 1);
        EXPECT_EQ(pos[0].second, 3);
        EXPECT_EQ(pos[1].second, 4);
        auto res = calcTrace(std::move(t1), pos);
        IContainer* tmp = res.get();
        EXPECT_NE(tmp, nullptr);
        EXPECT_NE(dynamic_cast<VectorContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<SparseContainer*>(tmp), nullptr);
        EXPECT_EQ(dynamic_cast<ScalarContainer*>(tmp), nullptr);
        EXPECT_EQ(res->rank(), 2);
        EXPECT_DOUBLE_EQ(res->element({4, 4}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({4, 4}).imag(), 0.1);
        EXPECT_DOUBLE_EQ(res->element({3, 2}).real(), 0.2);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).real(), 0.5);
    }

    TEST(TraceTest, Outer) {
        // Compare to direct evaluation
        auto amp = makeVector({2, 2, 11}, {SPIN_NUTAU, SPIN_TAUP, WILSON_BCTAUNU},
                              {(218526.018848 - 38094.3734443 * 1i),
                               (-41834.873825 + 0. * 1i),
                               0. * 1i,
                               (41834.873825 + 0. * 1i),
                               0. * 1i,
                               (-218526.018848 - 38094.3734443 * 1i),
                               0. * 1i,
                               (218526.018848 - 38094.3734443 * 1i),
                               0. * 1i,
                               (1994452.81665 + 625211.034433 * 1i),
                               0. * 1i,
                               (130542.449024 - 322278.778667 * 1i),
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               (69836.9021198 + 253350.864021 * 1i),
                               0. * 1i,
                               (130542.449024 - 322278.778667 * 1i),
                               0. * 1i,
                               (-432564.098293 - 609539.450478 * 1i),
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               (-130542.449024 - 322278.778667 * 1i),
                               0. * 1i,
                               (-69836.9021198 + 253350.864021 * 1i),
                               0. * 1i,
                               (432564.098293 - 609539.450478 * 1i),
                               0. * 1i,
                               0. * 1i,
                               (41834.873825 + 0. * 1i),
                               0. * 1i,
                               (-41834.873825 + 0. * 1i),
                               0. * 1i,
                               (218526.018848 + 38094.3734443 * 1i),
                               0. * 1i,
                               (-218526.018848 + 38094.3734443 * 1i),
                               0. * 1i,
                               (1994452.81665 - 625211.034433 * 1i)});
        EXPECT_DOUBLE_EQ(amp->element({0, 0, 1}).real(), (-41834.873825));
        EXPECT_DOUBLE_EQ(amp->element({0, 1, 0}).real(), (130542.449024));
        EXPECT_DOUBLE_EQ(amp->element({0, 1, 0}).imag(), (-322278.778667));
        EXPECT_DOUBLE_EQ(amp->element({1, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(amp->element({1, 0, 0}).imag(), 0.);
        auto t1 = makeOuterSquare(amp);
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 1);
        EXPECT_EQ(pos[1].first, 0);
        EXPECT_EQ(pos[0].second, 4);
        EXPECT_EQ(pos[1].second, 3);
        auto res = calcTrace(std::move(t1), pos);
        EXPECT_EQ(res->rank(), 2);
        auto amp2 = makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                               {1701.0974437794703,
                                -91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                -1188.3536642139165 - 389.30770234120007 * 1.i,
                                0.,
                                1701.0974437794703,
                                0.,
                                5519.964640738945 + 63.7469322621705 * 1.i,
                                0.,
                                -91.42008425985652 - 15.936733064847212 * 1.i,
                                17.501566679536705,
                                0.,
                                -17.501566679536705,
                                0.,
                                91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                -91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                -834.3768193446862 + 261.55624739502287 * 1.i,
                                0.,
                                0.,
                                0.,
                                17.501566679536705,
                                0.,
                                -17.501566679536705,
                                0.,
                                91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                -91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                834.3768193446862 + 261.55624739502287 * 1.i,
                                91.42008425985652 + 15.936733064847212 * 1.i,
                                -17.501566679536705,
                                0.,
                                17.501566679536705,
                                0.,
                                -91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                834.3768193446862 - 261.55624739502287 * 1.i,
                                0.,
                                0.,
                                0.,
                                -17.501566679536705,
                                0.,
                                17.501566679536705,
                                0.,
                                -91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                -834.3768193446862 - 261.55624739502287 * 1.i,
                                -1188.3536642139165 + 389.30770234120007 * 1.i,
                                91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                -91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                1182.6865539954797,
                                0.,
                                -1188.3536642139165 + 389.30770234120007 * 1.i,
                                0.,
                                -6442.931394564071 - 63.746932258253196 * 1.i,
                                0.,
                                0.,
                                0.,
                                91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                -91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                1701.0974437794703,
                                0.,
                                -1188.3536642139165 + 389.30770234120007 * 1.i,
                                0.,
                                5519.964640738945 - 63.7469322621705 * 1.i,
                                1701.0974437794703,
                                -91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                -1188.3536642139165 - 389.30770234120007 * 1.i,
                                0.,
                                1701.0974437794703,
                                0.,
                                5519.964640738945 + 63.7469322621705 * 1.i,
                                0.,
                                0.,
                                0.,
                                -91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                -1188.3536642139165 - 389.30770234120007 * 1.i,
                                0.,
                                1182.6865539954797,
                                0.,
                                -6442.931394564071 + 63.746932258253196 * 1.i,
                                5519.964640738945 - 63.74693226217096 * 1.i,
                                -834.3768193446862 - 261.5562473950229 * 1.i,
                                0.,
                                834.3768193446862 + 261.5562473950229 * 1.i,
                                0.,
                                -6442.93139456407 + 63.746932258253196 * 1.i,
                                0.,
                                5519.964640738945 - 63.74693226217096 * 1.i,
                                0.,
                                49273.80916240958,
                                0.,
                                0.,
                                0.,
                                834.3768193446862 - 261.5562473950229 * 1.i,
                                0.,
                                -834.3768193446862 + 261.5562473950229 * 1.i,
                                0.,
                                5519.964640738945 + 63.74693226217096 * 1.i,
                                0.,
                                -6442.93139456407 - 63.746932258253196 * 1.i,
                                0.,
                                49273.80916240958});
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_NEAR(const_cast<const IContainer*>(res.get())->element({idx1, idx2}).real() * 1.e-8,
                            amp2->element({idx1, idx2}).real(), 1.e-12);
                EXPECT_NEAR(const_cast<const IContainer*>(res.get())->element({idx1, idx2}).imag() * 1.e-8,
                            amp2->element({idx1, idx2}).imag(), 1.e-12);
            }
        }
    }

    TEST(TraceTest, Outer2) {
        // Compare to direct evaluation
        auto amp = makeVector({2, 2, 11}, {SPIN_NUTAU, SPIN_TAUP, WILSON_BCTAUNU},
                              {(218526.018848 - 38094.3734443 * 1i),
                               (-41834.873825 + 0. * 1i),
                               0. * 1i,
                               (41834.873825 + 0. * 1i),
                               0. * 1i,
                               (-218526.018848 - 38094.3734443 * 1i),
                               0. * 1i,
                               (218526.018848 - 38094.3734443 * 1i),
                               0. * 1i,
                               (1994452.81665 + 625211.034433 * 1i),
                               0. * 1i,
                               (130542.449024 - 322278.778667 * 1i),
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               (69836.9021198 + 253350.864021 * 1i),
                               0. * 1i,
                               (130542.449024 - 322278.778667 * 1i),
                               0. * 1i,
                               (-432564.098293 - 609539.450478 * 1i),
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               0. * 1i,
                               (-130542.449024 - 322278.778667 * 1i),
                               0. * 1i,
                               (-69836.9021198 + 253350.864021 * 1i),
                               0. * 1i,
                               (432564.098293 - 609539.450478 * 1i),
                               0. * 1i,
                               0. * 1i,
                               (41834.873825 + 0. * 1i),
                               0. * 1i,
                               (-41834.873825 + 0. * 1i),
                               0. * 1i,
                               (218526.018848 + 38094.3734443 * 1i),
                               0. * 1i,
                               (-218526.018848 + 38094.3734443 * 1i),
                               0. * 1i,
                               (1994452.81665 - 625211.034433 * 1i)});
        EXPECT_DOUBLE_EQ(amp->element({0, 0, 1}).real(), (-41834.873825));
        EXPECT_DOUBLE_EQ(amp->element({0, 1, 0}).real(), (130542.449024));
        EXPECT_DOUBLE_EQ(amp->element({0, 1, 0}).imag(), (-322278.778667));
        EXPECT_DOUBLE_EQ(amp->element({1, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(amp->element({1, 0, 0}).imag(), 0.);
        auto ampS = toSparse(std::move(amp));
        auto t1 = makeOuterSquare(ampS);
        auto pos = t1->getSpinLabelPairs();
        EXPECT_EQ(pos.size(), 2);
        EXPECT_EQ(pos[0].first, 1);
        EXPECT_EQ(pos[1].first, 0);
        EXPECT_EQ(pos[0].second, 4);
        EXPECT_EQ(pos[1].second, 3);
        auto res = calcTrace(std::move(t1), pos);
        EXPECT_EQ(res->rank(), 2);
        auto amp2 = makeVector({11, 11}, {WILSON_BCTAUNU, WILSON_BCTAUNU_HC},
                               {1701.0974437794703,
                                -91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                -1188.3536642139165 - 389.30770234120007 * 1.i,
                                0.,
                                1701.0974437794703,
                                0.,
                                5519.964640738945 + 63.7469322621705 * 1.i,
                                0.,
                                -91.42008425985652 - 15.936733064847212 * 1.i,
                                17.501566679536705,
                                0.,
                                -17.501566679536705,
                                0.,
                                91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                -91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                -834.3768193446862 + 261.55624739502287 * 1.i,
                                0.,
                                0.,
                                0.,
                                17.501566679536705,
                                0.,
                                -17.501566679536705,
                                0.,
                                91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                -91.42008425985652 - 15.936733064847212 * 1.i,
                                0.,
                                834.3768193446862 + 261.55624739502287 * 1.i,
                                91.42008425985652 + 15.936733064847212 * 1.i,
                                -17.501566679536705,
                                0.,
                                17.501566679536705,
                                0.,
                                -91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                834.3768193446862 - 261.55624739502287 * 1.i,
                                0.,
                                0.,
                                0.,
                                -17.501566679536705,
                                0.,
                                17.501566679536705,
                                0.,
                                -91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                91.42008425985652 + 15.936733064847212 * 1.i,
                                0.,
                                -834.3768193446862 - 261.55624739502287 * 1.i,
                                -1188.3536642139165 + 389.30770234120007 * 1.i,
                                91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                -91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                1182.6865539954797,
                                0.,
                                -1188.3536642139165 + 389.30770234120007 * 1.i,
                                0.,
                                -6442.931394564071 - 63.746932258253196 * 1.i,
                                0.,
                                0.,
                                0.,
                                91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                -91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                1701.0974437794703,
                                0.,
                                -1188.3536642139165 + 389.30770234120007 * 1.i,
                                0.,
                                5519.964640738945 - 63.7469322621705 * 1.i,
                                1701.0974437794703,
                                -91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                -1188.3536642139165 - 389.30770234120007 * 1.i,
                                0.,
                                1701.0974437794703,
                                0.,
                                5519.964640738945 + 63.7469322621705 * 1.i,
                                0.,
                                0.,
                                0.,
                                -91.42008425985654 + 15.936733064847214 * 1.i,
                                0.,
                                91.42008425985654 - 15.936733064847214 * 1.i,
                                0.,
                                -1188.3536642139165 - 389.30770234120007 * 1.i,
                                0.,
                                1182.6865539954797,
                                0.,
                                -6442.931394564071 + 63.746932258253196 * 1.i,
                                5519.964640738945 - 63.74693226217096 * 1.i,
                                -834.3768193446862 - 261.5562473950229 * 1.i,
                                0.,
                                834.3768193446862 + 261.5562473950229 * 1.i,
                                0.,
                                -6442.93139456407 + 63.746932258253196 * 1.i,
                                0.,
                                5519.964640738945 - 63.74693226217096 * 1.i,
                                0.,
                                49273.80916240958,
                                0.,
                                0.,
                                0.,
                                834.3768193446862 - 261.5562473950229 * 1.i,
                                0.,
                                -834.3768193446862 + 261.5562473950229 * 1.i,
                                0.,
                                5519.964640738945 + 63.74693226217096 * 1.i,
                                0.,
                                -6442.93139456407 - 63.746932258253196 * 1.i,
                                0.,
                                49273.80916240958});
        for (IndexType idx1 = 0; idx1 < 11; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 11; ++idx2) {
                EXPECT_NEAR(const_cast<const IContainer*>(res.get())->element({idx1, idx2}).real() * 1.e-8,
                            amp2->element({idx1, idx2}).real(), 1.e-12);
                EXPECT_NEAR(const_cast<const IContainer*>(res.get())->element({idx1, idx2}).imag() * 1.e-8,
                            amp2->element({idx1, idx2}).imag(), 1.e-12);
            }
        }
    }
} // namespace Hammer::MultiDimensional
