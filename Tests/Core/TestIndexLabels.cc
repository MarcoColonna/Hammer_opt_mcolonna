///
/// @file  TestIndexLabels.cc
/// @brief Tests for IndexLabel classification and manipulation functions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#include "Hammer/IndexLabels.hh"

namespace Hammer {

    TEST(IndexLabelsTest, IsSpinIndexPositive) {
        EXPECT_TRUE(isSpinIndex(SPIN_TAUP));
        EXPECT_TRUE(isSpinIndex(SPIN_DSTAR));
        EXPECT_TRUE(isSpinIndex(SPIN_JPSI));
        EXPECT_TRUE(isSpinIndex(SPIN_LCSTAR32));
    }

    TEST(IndexLabelsTest, IsSpinIndexNegativeHC) {
        EXPECT_TRUE(isSpinIndex(SPIN_TAUP_HC));
        EXPECT_TRUE(isSpinIndex(SPIN_DSTAR_HC));
        EXPECT_TRUE(isSpinIndex(SPIN_LCSTAR32_HC));
    }

    TEST(IndexLabelsTest, IsSpinIndexNone) {
        EXPECT_TRUE(isSpinIndex(NONE));
    }

    TEST(IndexLabelsTest, IsSpinIndexBoundary) {
        EXPECT_TRUE(isSpinIndex(static_cast<IndexLabel>(SPIN_INDEX_END)));
        EXPECT_FALSE(isSpinIndex(static_cast<IndexLabel>(SPIN_INDEX_END + 1)));
    }

    TEST(IndexLabelsTest, IsSpinIndexRefLabels) {
        EXPECT_TRUE(isSpinIndex(SPIN_TAUP_REF));
        EXPECT_TRUE(isSpinIndex(SPIN_NUMU_BAR_REF));
        EXPECT_TRUE(isSpinIndex(SPIN_RHO));
        EXPECT_TRUE(isSpinIndex(SPIN_LB));
    }

    TEST(IndexLabelsTest, IsFFIndexPositive) {
        EXPECT_TRUE(isFFIndex(FF_BD));
        EXPECT_TRUE(isFFIndex(FF_BDSTAR));
        EXPECT_TRUE(isFFIndex(FF_BCJPSI));
        EXPECT_TRUE(isFFIndex(FF_LBLC));
        EXPECT_TRUE(isFFIndex(FF_BPI));
        EXPECT_TRUE(isFFIndex(FF_BRHO));
    }

    TEST(IndexLabelsTest, IsFFIndexVarRange) {
        EXPECT_TRUE(isFFIndex(FF_BD_VAR));
        EXPECT_TRUE(isFFIndex(FF_BDSTAR_VAR));
        EXPECT_TRUE(isFFIndex(FF_BCJPSI_VAR));
        EXPECT_TRUE(isFFIndex(FF_LBLC_VAR));
    }

    TEST(IndexLabelsTest, IsFFIndexNegativeHC) {
        EXPECT_TRUE(isFFIndex(FF_BD_HC));
        EXPECT_TRUE(isFFIndex(FF_BDSTAR_HC));
        EXPECT_TRUE(isFFIndex(FF_BD_VAR_HC));
    }

    TEST(IndexLabelsTest, IsFFIndexBoundaries) {
        EXPECT_FALSE(isFFIndex(static_cast<IndexLabel>(FF_INDEX_START - 1)));
        EXPECT_TRUE(isFFIndex(static_cast<IndexLabel>(FF_INDEX_START)));
        EXPECT_TRUE(isFFIndex(static_cast<IndexLabel>(FF_VAR_INDEX_END)));
        EXPECT_FALSE(isFFIndex(static_cast<IndexLabel>(FF_VAR_INDEX_END + 1)));
    }

    TEST(IndexLabelsTest, IsFFIndexExcludesSpinAndWC) {
        EXPECT_FALSE(isFFIndex(SPIN_TAUP));
        EXPECT_FALSE(isFFIndex(NONE));
        EXPECT_FALSE(isFFIndex(WILSON_BCTAUNU));
    }

    TEST(IndexLabelsTest, IsWCIndexPositive) {
        EXPECT_TRUE(isWCIndex(WILSON_BCTAUNU));
        EXPECT_TRUE(isWCIndex(WILSON_BCMUNU));
        EXPECT_TRUE(isWCIndex(WILSON_BCENU));
        EXPECT_TRUE(isWCIndex(WILSON_BUTAUNU));
        EXPECT_TRUE(isWCIndex(WILSON_BUMUNU));
        EXPECT_TRUE(isWCIndex(WILSON_BUENU));
    }

    TEST(IndexLabelsTest, IsWCIndexNegativeHC) {
        EXPECT_TRUE(isWCIndex(WILSON_BCTAUNU_HC));
        EXPECT_TRUE(isWCIndex(WILSON_BCMUNU_HC));
    }

    TEST(IndexLabelsTest, IsWCIndexBoundaries) {
        EXPECT_FALSE(isWCIndex(static_cast<IndexLabel>(WC_INDEX_START - 1)));
        EXPECT_TRUE(isWCIndex(static_cast<IndexLabel>(WC_INDEX_START)));
        EXPECT_TRUE(isWCIndex(static_cast<IndexLabel>(WC_INDEX_END)));
        EXPECT_FALSE(isWCIndex(static_cast<IndexLabel>(WC_INDEX_END + 1)));
    }

    TEST(IndexLabelsTest, IsWCIndexExcludesOthers) {
        EXPECT_FALSE(isWCIndex(NONE));
        EXPECT_FALSE(isWCIndex(SPIN_TAUP));
        EXPECT_FALSE(isWCIndex(FF_BD));
    }

    TEST(IndexLabelsTest, IsSpecializedIndexFalseForPlainLabels) {
        EXPECT_FALSE(isSpecializedIndex(WILSON_BCTAUNU));
        EXPECT_FALSE(isSpecializedIndex(SPIN_TAUP));
        EXPECT_FALSE(isSpecializedIndex(FF_BD));
        EXPECT_FALSE(isSpecializedIndex(NONE));
    }

    TEST(IndexLabelsTest, IsSpecializedIndexTrueForLargeLabel) {
        auto specLabel = specializeLabel(WILSON_BCTAUNU, 1u);
        EXPECT_TRUE(isSpecializedIndex(specLabel));
    }

    TEST(IndexLabelsTest, SpecializeLabelPositiveIndex) {
        IndexLabel base = WILSON_BCTAUNU;
        auto specLabel = specializeLabel(base, 1u);
        EXPECT_GT(static_cast<int64_t>(specLabel), static_cast<int64_t>(SPECIALIZATION_PAD));
        EXPECT_TRUE(isSpecializedIndex(specLabel));
    }

    TEST(IndexLabelsTest, SpecializeLabelNegativeIndex) {
        IndexLabel base = WILSON_BCTAUNU_HC;
        auto specLabel = specializeLabel(base, 1u);
        EXPECT_LT(static_cast<int64_t>(specLabel), -static_cast<int64_t>(SPECIALIZATION_PAD));
        EXPECT_TRUE(isSpecializedIndex(specLabel));
    }

    TEST(IndexLabelsTest, GetSpecializedLabelComponents) {
        uint32_t pad = 2u;
        IndexLabel base = WILSON_BCTAUNU;
        auto specLabel = specializeLabel(base, pad);
        auto [recoveredPad, recoveredBase] = getSpecializedLabelComponents(specLabel);
        EXPECT_EQ(recoveredPad, pad);
        EXPECT_EQ(recoveredBase, static_cast<int>(base));
    }

    TEST(IndexLabelsTest, GetSpecializedLabelComponentsNegative) {
        uint32_t pad = 3u;
        IndexLabel base = WILSON_BCTAUNU_HC;
        auto specLabel = specializeLabel(base, pad);
        auto [recoveredPad, recoveredBase] = getSpecializedLabelComponents(specLabel);
        EXPECT_EQ(recoveredPad, pad);
        EXPECT_EQ(recoveredBase, static_cast<int>(base));
    }

    TEST(IndexLabelsTest, GenerateSpecializedPadDifferentStrings) {
        auto pad1 = generateSpecializedPad("SpecA");
        auto pad2 = generateSpecializedPad("SpecB");
        EXPECT_NE(pad1, pad2);
    }

    TEST(IndexLabelsTest, GenerateSpecializedPadDeterministic) {
        auto pad1 = generateSpecializedPad("TestSpec");
        auto pad2 = generateSpecializedPad("TestSpec");
        EXPECT_EQ(pad1, pad2);
    }

    TEST(IndexLabelsTest, SpecializeLabelRoundTripWithPad) {
        uint32_t pad = generateSpecializedPad("RH");
        IndexLabel base = WILSON_BCMUNU;
        auto specLabel = specializeLabel(base, pad);

        EXPECT_TRUE(isSpecializedIndex(specLabel));
        EXPECT_FALSE(isSpinIndex(specLabel));
        EXPECT_FALSE(isFFIndex(specLabel));
        EXPECT_FALSE(isWCIndex(specLabel));

        auto [recoveredPad, recoveredBase] = getSpecializedLabelComponents(specLabel);
        EXPECT_EQ(recoveredPad, pad);
        EXPECT_EQ(recoveredBase, static_cast<int>(base));
    }

    TEST(IndexLabelsTest, IntegrationIndexNotInAnyStandardRange) {
        auto label = static_cast<IndexLabel>(INTEGRATION_INDEX);
        EXPECT_FALSE(isSpinIndex(label));
        EXPECT_FALSE(isWCIndex(label));
        EXPECT_FALSE(isFFIndex(label));
        EXPECT_FALSE(isSpecializedIndex(label));
    }

} // namespace Hammer
