///
/// @file  TestBlockIndexing.cc
/// @brief Tests for BlockIndexing
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/BlockIndexing.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    class TestableBlockIndexing : public BlockIndexing {
    public:

        TestableBlockIndexing() = default;
        TestableBlockIndexing(const std::vector<IndexList>& dims, const std::vector<LabelsList>& labels)
            : BlockIndexing{dims, labels} {
        }
        TestableBlockIndexing(const LabeledIndexing<AlignedIndexing>& left,
                              const LabeledIndexing<AlignedIndexing>& right)
            : BlockIndexing{left, right} {
        }
        ~TestableBlockIndexing() = default;
    };

    TEST(BlockTest, ConstructionAccess1) {
        vector<IndexList> v1{{2, 3}, {4, 1, 7}, {11}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_DSSD2STAR, FF_BD, WILSON_BCENU}, {INTEGRATION_INDEX}};
        TestableBlockIndexing t{v1, v2};
        EXPECT_EQ(t.rank(), 6);
        EXPECT_EQ(t.dim(0), 2);
        EXPECT_EQ(t.dim(1), 3);
        EXPECT_EQ(t.dim(2), 4);
        EXPECT_EQ(t.dim(3), 1);
        EXPECT_EQ(t.dim(4), 7);
        EXPECT_EQ(t.dim(5), 11);
        EXPECT_EQ(t.dims()[3], 1);
        EXPECT_EQ(t.labels()[3], FF_BD);
        EXPECT_EQ(t.dim(SPIN_TAUM), 2);
        EXPECT_EQ(t.dim(SPIN_DSTAR), 3);
        EXPECT_EQ(t.dim(SPIN_DSSD2STAR), 4);
        EXPECT_EQ(t.dim(FF_BD), 1);
        EXPECT_EQ(t.dim(WILSON_BCENU), 7);
        EXPECT_EQ(t.dim(INTEGRATION_INDEX), 11);
        EXPECT_EQ(t.numValues(), 2 * 3 * 4 * 7 * 11);
        EXPECT_EQ(t.spinIndices().size(), 3);
    }

    TEST(BlockTest, ConstructionAccess2) {
        LabeledIndexing<AlignedIndexing> id1{{2, 2, 11}, {SPIN_DSTAR, SPIN_TAUM, WILSON_BCENU}};
        LabeledIndexing<AlignedIndexing> id2{{2, 2, 11}, {SPIN_DSTAR_HC, SPIN_TAUM_HC, WILSON_BCENU_HC}};
        TestableBlockIndexing t{id1, id2};
        EXPECT_EQ(t.rank(), 6);
        EXPECT_EQ(t.dim(0), 2);
        EXPECT_EQ(t.dim(1), 2);
        EXPECT_EQ(t.dim(2), 11);
        EXPECT_EQ(t.dim(3), 2);
        EXPECT_EQ(t.dim(4), 2);
        EXPECT_EQ(t.dim(5), 11);
        EXPECT_EQ(t.dims()[3], 2);
        EXPECT_EQ(t.labels()[3], SPIN_DSTAR_HC);
        EXPECT_EQ(t.dim(SPIN_DSTAR), 2);
        EXPECT_EQ(t.dim(SPIN_TAUM), 2);
        EXPECT_EQ(t.dim(WILSON_BCENU), 11);
        EXPECT_THROW((void) t.dim(FF_BD), IndexLabelError);
        EXPECT_EQ(t.dim(WILSON_BCENU), 11);
        EXPECT_EQ(t.numValues(), 4 * 4 * 121);
        EXPECT_EQ(t.spinIndices().size(), 4);
    }

    TEST(BlockTest, LabelAccess) {
        vector<IndexList> v1{{2, 3}, {4, 1, 7}, {11}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_DSSD2STAR, FF_BD, WILSON_BCENU}, {INTEGRATION_INDEX}};
        TestableBlockIndexing t{v1, v2};
        EXPECT_EQ(t.labelIndex(SPIN_TAUM), 0);
        EXPECT_EQ(t.labelIndex(WILSON_BCENU), 4);
        EXPECT_EQ(t.labelIndex(FF_BD_VAR), 6);
    }

    TEST(BlockTest, IndexAccess) {
        vector<IndexList> v1{{2, 3}, {4, 1, 7}, {11}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_DSSD2STAR, FF_BD, WILSON_BCENU}, {INTEGRATION_INDEX}};
        TestableBlockIndexing t{v1, v2};
        EXPECT_EQ(t.getElementIndex(3).first, 1);
        EXPECT_EQ(t.getElementIndex(3).second, 1);
        EXPECT_EQ(t.getElementIndex(5).first, 2);
        EXPECT_EQ(t.getElementIndex(5).second, 0);
    }

    TEST(BlockTest, ValidIndices) {
        vector<IndexList> v1{{2, 3}, {4, 1, 7}, {11}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_DSSD2STAR, FF_BD, WILSON_BCENU}, {INTEGRATION_INDEX}};
        TestableBlockIndexing t{v1, v2};
        EXPECT_TRUE(t.checkValidIndices({1, 2, 0, 0, 6, 5}));
        EXPECT_FALSE(t.checkValidIndices({1, 2, 0, 0, 6, 15}));
        EXPECT_TRUE(t.checkValidIndices({{1, 2}, {0, 0, 6}, {5}}));
        EXPECT_FALSE(t.checkValidIndices({{1, 2, 0}, {0, 6}, {5}}));
    }

    TEST(BlockTest, SameLabels) {
        vector<IndexList> v1{{2, 3}, {4, 1, 7}, {11}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_DSSD2STAR, FF_BD, WILSON_BCENU}, {INTEGRATION_INDEX}};
        TestableBlockIndexing t{v1, v2};
        auto res = t.getSameLabelPairs({FF_BD, SPIN_TAUM, SPIN_DSSD1, WILSON_BCTAUNU},
                                       {SPIN_TAUM, SPIN_DSSD1, FF_BD, WILSON_BCTAUNU});
        EXPECT_EQ(res.size(), 2);
        EXPECT_EQ(res[0].first, 3);
        EXPECT_EQ(res[0].second, 0);
        EXPECT_EQ(res[1].first, 0);
        EXPECT_EQ(res[1].second, 1);
    }

    TEST(BlockTest, OppositeLabels) {
        vector<IndexList> v1{{2, 3}, {2, 3, 7}, {11}};
        vector<LabelsList> v2{
            {SPIN_TAUM, SPIN_DSTAR}, {SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCENU}, {INTEGRATION_INDEX}};
        TestableBlockIndexing t{v1, v2};
        auto res = t.getOppositeLabelPairs({SPIN_TAUM, SPIN_DSTAR, WILSON_BCENU, INTEGRATION_INDEX});
        EXPECT_EQ(res.size(), 2);
        EXPECT_EQ(res[0].first, 0);
        EXPECT_EQ(res[0].second, 2);
        EXPECT_EQ(res[1].first, 1);
        EXPECT_EQ(res[1].second, 3);
    }

    TEST(BlockTest, SameShape) {
        vector<IndexList> v1{{2, 3}, {4, 1, 7}, {11}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_DSSD2STAR, FF_BD, WILSON_BCENU}, {INTEGRATION_INDEX}};
        vector<IndexList> v3{{2, 3, 4}, {1, 7}, {11}};
        vector<LabelsList> v4{{SPIN_TAUM, SPIN_DSTAR, SPIN_DSSD2STAR}, {FF_BD, WILSON_BCENU}, {INTEGRATION_INDEX}};
        vector<IndexList> v5{{2, 3, 8}, {1, 7}, {11}};
        TestableBlockIndexing t{v1, v2};
        TestableBlockIndexing t2{v1, v2};
        TestableBlockIndexing t3{v3, v4};
        TestableBlockIndexing t4{v5, v4};
        LabeledIndexing<AlignedIndexing> id1{
            {2, 3, 4, 1, 7, 11}, {SPIN_TAUM, SPIN_DSTAR, SPIN_DSSD2STAR, FF_BD, WILSON_BCENU, INTEGRATION_INDEX}};
        LabeledIndexing<AlignedIndexing> id2{
            {2, 3, 4, 8, 7, 11}, {SPIN_TAUM, SPIN_DSTAR, SPIN_DSSD2STAR, FF_BD, WILSON_BCENU, INTEGRATION_INDEX}};
        EXPECT_TRUE(t.isSameLabelShape(t2));
        EXPECT_TRUE(t.isSameLabelShape(t3));
        EXPECT_FALSE(t.isSameLabelShape(t3, true));
        EXPECT_FALSE(t.isSameLabelShape(t4));
        EXPECT_TRUE(t.isSameLabelShape(id1));
        EXPECT_FALSE(t.isSameLabelShape(id2));
    }

    TEST(BlockTest, Split) {
        vector<IndexList> v1{{2, 3}, {2, 3, 7}, {11}};
        vector<LabelsList> v2{
            {SPIN_TAUM, SPIN_DSTAR}, {SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCENU}, {INTEGRATION_INDEX}};
        TestableBlockIndexing t{v1, v2};
        auto res = t.splitIndices({1, 2, 0, 2, 5, 8});
        EXPECT_EQ(res.size(), 3);
        EXPECT_EQ(res[0][0], 1);
        EXPECT_EQ(res[0][1], 2);
        EXPECT_EQ(res[1][0], 0);
        EXPECT_EQ(res[1][1], 2);
        EXPECT_EQ(res[1][2], 5);
        EXPECT_EQ(res[2][0], 8);
    }

    TEST(BlockTest, ProcessShifts) {
        vector<IndexList> v1{{2, 3}, {2, 3, 7}, {9}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCMUNU}, {WILSON_BCENU}};
        vector<IndexList> v3{{2, 1}, {5, 3, 9}, {11, 5}};
        vector<LabelsList> v4{
            {SPIN_TAUM_HC, SPIN_DSSD1}, {SPIN_DSTAR_HC, SPIN_DSTAR, WILSON_BCENU}, {INTEGRATION_INDEX, WILSON_BCTAUNU}};
        TestableBlockIndexing t1{v1, v2};
        TestableBlockIndexing t2{v3, v4};
        DotGroupList chks{DotGroupType{{}, {2}, {}}, DotGroupType{{1}, {0}, {{2, 0}}},
                          DotGroupType{{0, 2}, {1}, {{1, 3}, {5, 4}}}};
        auto res = t1.processShifts(chks, IndexPairMember::Left);
        EXPECT_EQ(res.size(), 2);
        EXPECT_EQ(get<0>(res[0]).size(), 6);
        EXPECT_EQ(get<0>(res[0])[5], 12);
        EXPECT_EQ(get<0>(res[0])[4], 0);
        EXPECT_EQ(get<0>(res[0])[3], 3);
        EXPECT_EQ(get<0>(res[0])[2], 0);
        EXPECT_EQ(get<0>(res[0])[1], 12);
        EXPECT_EQ(get<0>(res[0])[0], 12);
        EXPECT_TRUE(get<1>(res[0])[5]);
        EXPECT_TRUE(get<1>(res[0])[4]);
        EXPECT_TRUE(get<1>(res[0])[3]);
        EXPECT_FALSE(get<1>(res[0])[2]);
        EXPECT_TRUE(get<1>(res[0])[1]);
        EXPECT_TRUE(get<1>(res[0])[0]);
        EXPECT_EQ(get<2>(res[0]), 4096 / 2 / 4 / 16);
        // second element
        EXPECT_EQ(get<0>(res[1]).size(), 6);
        EXPECT_EQ(get<0>(res[1]).size(), 6);
        EXPECT_EQ(get<0>(res[1])[5], 1);
        EXPECT_EQ(get<0>(res[1])[4], 7);
        EXPECT_EQ(get<0>(res[1])[3], 7);
        EXPECT_EQ(get<0>(res[1])[2], 7);
        EXPECT_EQ(get<0>(res[1])[1], 0);
        EXPECT_EQ(get<0>(res[1])[0], 0);
        EXPECT_FALSE(get<1>(res[1])[5]);
        EXPECT_TRUE(get<1>(res[1])[4]);
        EXPECT_TRUE(get<1>(res[1])[3]);
        EXPECT_TRUE(get<1>(res[1])[2]);
        EXPECT_FALSE(get<1>(res[1])[1]);
        EXPECT_TRUE(get<1>(res[1])[0]);
        EXPECT_EQ(get<2>(res[1]), 128 / 2 / 4 / 8);
    }

    TEST(BlockTest, ProcessShifts2) {
        vector<IndexList> v1{{3}, {3}};
        vector<LabelsList> v2{{WILSON_BCMUNU}, {WILSON_BCMUNU_HC}};
        vector<IndexList> v3{{2, 3}, {2, 3}};
        vector<LabelsList> v4{{FF_BD, WILSON_BCMUNU}, {FF_BD_HC, WILSON_BCMUNU_HC}};
        TestableBlockIndexing t1{v1, v2};
        TestableBlockIndexing t2{v3, v4};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        auto res1 = t1.processShifts(chks, IndexPairMember::Left);
        EXPECT_EQ(res1.size(), 2);
        // first element
        EXPECT_EQ(get<0>(res1[0]).size(), 2);
        EXPECT_EQ(get<0>(res1[0])[0], 0);
        EXPECT_EQ(get<0>(res1[0])[1], 2);
        EXPECT_FALSE(get<1>(res1[0])[0]);
        EXPECT_TRUE(get<1>(res1[0])[1]);
        EXPECT_EQ(get<2>(res1[0]), 4 / 4);
        // second element
        EXPECT_EQ(get<0>(res1[1]).size(), 2);
        EXPECT_EQ(get<0>(res1[1])[0], 2);
        EXPECT_EQ(get<0>(res1[1])[1], 0);
        EXPECT_TRUE(get<1>(res1[1])[0]);
        EXPECT_FALSE(get<1>(res1[1])[1]);
        EXPECT_EQ(get<2>(res1[1]), 4 / 4);
        auto res2 = t2.processShifts(chks, IndexPairMember::Right);
        EXPECT_EQ(res2.size(), 2);
        // first element
        EXPECT_EQ(get<0>(res2[0]).size(), 4);
        EXPECT_EQ(get<0>(res2[0])[0], 0);
        EXPECT_EQ(get<0>(res2[0])[1], 0);
        EXPECT_EQ(get<0>(res2[0])[2], 1 + 2 + 1);
        EXPECT_EQ(get<0>(res2[0])[3], 1 + 2 + 1);
        EXPECT_TRUE(get<1>(res2[0])[0]);
        EXPECT_FALSE(get<1>(res2[0])[1]);
        EXPECT_TRUE(get<1>(res2[0])[2]);
        EXPECT_TRUE(get<1>(res2[0])[3]);
        EXPECT_EQ(get<2>(res2[0]), 2 * 4 * 2 / (2 * 4));
        // second element
        EXPECT_EQ(get<0>(res2[1]).size(), 4);
        EXPECT_EQ(get<0>(res2[1])[0], 1 + 2 + 1); // sum of pads of spectator indices in the full tensor
        EXPECT_EQ(get<0>(res2[1])[1], 1 + 2 + 1);
        EXPECT_EQ(get<0>(res2[1])[2], 0);
        EXPECT_EQ(get<0>(res2[1])[3], 0);
        EXPECT_TRUE(get<1>(res2[1])[0]);
        EXPECT_TRUE(get<1>(res2[1])[1]);
        EXPECT_TRUE(get<1>(res2[1])[2]);
        EXPECT_FALSE(get<1>(res2[1])[3]);
        EXPECT_EQ(get<2>(res2[1]), 2 * 4 * 2 / (2 * 4));
    }

    TEST(BlockTest, ProcessShifts3) {
        vector<IndexList> v1{{2}, {2}};
        vector<LabelsList> v2{{WILSON_BCMUNU}, {WILSON_BCMUNU_HC}};
        vector<IndexList> v3{{3, 2}, {3, 2}};
        vector<LabelsList> v4{{FF_BD, WILSON_BCMUNU}, {FF_BD_HC, WILSON_BCMUNU_HC}};
        TestableBlockIndexing t1{v1, v2};
        TestableBlockIndexing t2{v3, v4};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        auto res1 = t1.processShifts(chks, IndexPairMember::Left);
        EXPECT_EQ(res1.size(), 2);
        // first element
        EXPECT_EQ(get<0>(res1[0]).size(), 2);
        EXPECT_EQ(get<0>(res1[0])[0], 0);
        EXPECT_EQ(get<0>(res1[0])[1], 1);
        EXPECT_FALSE(get<1>(res1[0])[0]);
        EXPECT_TRUE(get<1>(res1[0])[1]);
        EXPECT_EQ(get<2>(res1[0]), 2 / 2);
        // second element
        EXPECT_EQ(get<0>(res1[1]).size(), 2);
        EXPECT_EQ(get<0>(res1[1])[0], 1);
        EXPECT_EQ(get<0>(res1[1])[1], 0);
        EXPECT_TRUE(get<1>(res1[1])[0]);
        EXPECT_FALSE(get<1>(res1[1])[1]);
        EXPECT_EQ(get<2>(res1[1]), 2 / 2);
        auto res2 = t2.processShifts(chks, IndexPairMember::Right);
        EXPECT_EQ(res2.size(), 2);
        // first element
        EXPECT_EQ(get<0>(res2[0]).size(), 4);
        EXPECT_EQ(get<0>(res2[0])[0], 0);
        EXPECT_EQ(get<0>(res2[0])[1], 0);
        EXPECT_EQ(get<0>(res2[0])[2], 1 + 2 + 2);
        EXPECT_EQ(get<0>(res2[0])[3], 1 + 2 + 2);
        EXPECT_TRUE(get<1>(res2[0])[0]);
        EXPECT_FALSE(get<1>(res2[0])[1]);
        EXPECT_TRUE(get<1>(res2[0])[2]);
        EXPECT_TRUE(get<1>(res2[0])[3]);
        EXPECT_EQ(get<2>(res2[0]), 2 * 4 * 4 / (2 * 4));
        // second element
        EXPECT_EQ(get<0>(res2[1]).size(), 4);
        EXPECT_EQ(get<0>(res2[1])[0], 1 + 2 + 2); // sum of pads of spectator indices in the full tensor
        EXPECT_EQ(get<0>(res2[1])[1], 1 + 2 + 2);
        EXPECT_EQ(get<0>(res2[1])[2], 0);
        EXPECT_EQ(get<0>(res2[1])[3], 0);
        EXPECT_TRUE(get<1>(res2[1])[0]);
        EXPECT_TRUE(get<1>(res2[1])[1]);
        EXPECT_TRUE(get<1>(res2[1])[2]);
        EXPECT_FALSE(get<1>(res2[1])[3]);
        EXPECT_EQ(get<2>(res2[1]), 2 * 4 * 4 / (2 * 4));
    }

    TEST(BlockTest, OuterIter) {
        vector<IndexList> v1{{2, 3}, {2, 3, 7}, {9}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCMUNU}, {WILSON_BCENU}};
        auto t1 = SharedTensorData{makeEmptySparse({2, 3}, {SPIN_TAUM, SPIN_DSTAR}).release()};
        auto t2 = SharedTensorData{makeEmptySparse({2, 3, 7}, {SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCMUNU}).release()};
        auto t3 = SharedTensorData{makeEmptySparse({9}, {WILSON_BCENU}).release()};
        OuterElemIterator::EntryType entry{{t1, false}, {t2, false}, {t3, true}};
        OuterElemIterator it{entry};
        auto itb = it.begin();
        auto ite = it.end();
        EXPECT_TRUE(it == itb);
        EXPECT_TRUE(itb == ite);
        t1->element({0, 1}) = 1.0;
        t1->element({1, 2}) = 2.0;
        t2->element({0, 2, 4}) = 3.0;
        t2->element({0, 2, 6}) = 4.0;
        t3->element({1}) = 0.1;
        t3->element({4}) = 0.2i;
        t3->element({7}) = 0.3i;
        t3->element({8}) = 0.4;
        OuterElemIterator it2{entry};
        auto itb2 = it2.begin();
        auto ite2 = it2.end();
        EXPECT_EQ(itb2.position(0), 1);
        EXPECT_EQ(itb2.position(1), (2 << 3) + 4);
        EXPECT_EQ(itb2.position(2), 1);
        ++it2;
        EXPECT_EQ(it2.position(0), 1);
        EXPECT_EQ(it2.position(1), (2 << 3) + 4);
        EXPECT_EQ(it2.position(2), 4);
        EXPECT_DOUBLE_EQ((*it2).real(), 0.);
        EXPECT_DOUBLE_EQ((*it2).imag(), -0.2 * 1.0 * 3.0);
        ++it2;
        EXPECT_EQ(it2.position(2), 7);
        ++it2;
        EXPECT_EQ(it2.position(2), 8);
        ++it2;
        EXPECT_EQ(it2.position(0), 1);
        EXPECT_EQ(it2.position(1), (2 << 3) + 6);
        EXPECT_EQ(it2.position(2), 1);
        ++it2;
        ++it2;
        ++it2;
        ++it2;
        EXPECT_EQ(it2.position(0), (1 << 2) + 2);
        EXPECT_EQ(it2.position(1), (2 << 3) + 4);
        EXPECT_EQ(it2.position(2), 1);
    }

    TEST(BlockTest, OuterIter2) {
        auto t1 = SharedTensorData{makeEmptySparse({2}, {SPIN_TAUM}).release()};
        OuterElemIterator::EntryType entry{{t1, false}};
        OuterElemIterator it{entry};
        auto itb = it.begin();
        auto ite = it.end();
        EXPECT_TRUE(it == itb);
        EXPECT_TRUE(itb == ite);
        t1->element({0}) = 1.0;
        t1->element({1}) = 2.0;
        OuterElemIterator it2{entry};
        auto itb2 = it2.begin();
        auto ite2 = it2.end();
        EXPECT_EQ(itb2.position(0), 0);
        ++it2;
        EXPECT_EQ(it2.position(0), 1);
        EXPECT_DOUBLE_EQ((*it2).real(), 2.);
        EXPECT_DOUBLE_EQ((*it2).imag(), 0.);
    }

    TEST(BlockTest, OuterIter3) {
        auto t1 = SharedTensorData{makeEmptySparse({3}, {SPIN_TAUM}).release()};
        OuterElemIterator::EntryType entry{{t1, false}, {t1, true}};
        OuterElemIterator it{entry};
        auto itb = it.begin();
        auto ite = it.end();
        EXPECT_TRUE(it == itb);
        EXPECT_TRUE(itb == ite);
        t1->element({0}) = 1.0;
        t1->element({1}) = 2.0;
        OuterElemIterator it2{entry};
        auto itb2 = it2.begin();
        auto ite2 = it2.end();
        EXPECT_EQ(itb2.position(0), 0);
        EXPECT_EQ(itb2.position(1), 0);
        ++it2;
        EXPECT_DOUBLE_EQ((*it2).real(), 2.);
        EXPECT_DOUBLE_EQ((*it2).imag(), 0.);
        EXPECT_EQ(it2.position(0), 0);
        EXPECT_EQ(it2.position(1), 1);
        ++it2;
        EXPECT_EQ(it2.position(0), 1);
        EXPECT_EQ(it2.position(1), 0);
        ++it2;
        EXPECT_EQ(it2.position(0), 1);
        EXPECT_EQ(it2.position(1), 1);
        ++it2;
        EXPECT_TRUE(it2 == ite2);
    }

    TEST(BlockTest, OuterIter4) {
        auto t1 = SharedTensorData{makeVector({5}, {WILSON_BCENU}, {1., 0., 0., 0., 0.}).release()};
        auto t2 = SharedTensorData{makeVector({7}, {FF_BD_VAR}, {1., 2., 3., 4., 5., 6., 7.}).release()};
        OuterElemIterator::EntryType entry{{t1, false}, {t1, true}, {t2, false}, {t2, true}};
        OuterElemIterator it{entry};
        auto itb = it.begin();
        auto ite = it.end();
        EXPECT_EQ(itb.position(0), 0);
        EXPECT_EQ(itb.position(1), 0);
        EXPECT_EQ(itb.position(2), 0);
        EXPECT_EQ(itb.position(3), 0);
        for (size_t i = 0; i < 7; ++i) {
            EXPECT_EQ(it.position(0), 0);
            EXPECT_EQ(it.position(1), 0);
            EXPECT_EQ(it.position(2), 0);
            EXPECT_EQ(it.position(3), i);
            ++it;
        }
        EXPECT_EQ(it.position(0), 0);
        EXPECT_EQ(it.position(1), 0);
        EXPECT_EQ(it.position(2), 1);
        EXPECT_EQ(it.position(3), 0);
    }

    TEST(BlockTest, FullPosition) {
        vector<IndexList> v1{{2, 3}, {2, 3, 7}, {9}};
        vector<LabelsList> v2{{SPIN_TAUM, SPIN_DSTAR}, {SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCMUNU}, {WILSON_BCENU}};
        auto t1 = SharedTensorData{makeEmptySparse({2, 3}, {SPIN_TAUM, SPIN_DSTAR}).release()};
        auto t2 = SharedTensorData{makeEmptySparse({2, 3, 7}, {SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCMUNU}).release()};
        auto t3 = SharedTensorData{makeEmptySparse({9}, {WILSON_BCENU}).release()};
        BlockIndexing bi{v1, v2};
        DotGroupList chks{DotGroupType{{}, {2}, {}}, DotGroupType{{1}, {0}, {{2, 0}}},
                          DotGroupType{{0, 2}, {1}, {{1, 3}, {5, 4}}}};
        t1->element({0, 1}) = 1.0;
        t1->element({1, 2}) = 2.0;
        t2->element({0, 2, 4}) = 3.0;
        t2->element({0, 2, 6}) = 4.0;
        t3->element({1}) = 0.1;
        t3->element({4}) = 0.2i;
        t3->element({7}) = 0.3i;
        t3->element({8}) = 0.4;
        OuterElemIterator::EntryType entry0{{t2, false}};
        OuterElemIterator it0{entry0};
        auto ite0 = it0.end();
        auto chk0 = get<0>(chks[1]);
        OuterElemIterator::EntryType entry1{{t1, false}, {t3, false}};
        OuterElemIterator it1{entry1};
        auto ite1 = it1.end();
        auto chk1 = get<0>(chks[2]);
        LabeledIndexing<AlignedIndexing> other{
            {2, 3, 2, 3, 7, 9}, {SPIN_TAUM, SPIN_DSTAR, SPIN_TAUM_HC, SPIN_DSTAR_HC, WILSON_BCMUNU, WILSON_BCENU}};
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 0, 0, 2, 4, 0}));
        ++it0;
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 0, 0, 2, 6, 0}));
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 1, 0, 0, 0, 1}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 1, 0, 0, 0, 4}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 1, 0, 0, 0, 7}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 1, 0, 0, 0, 8}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({1, 2, 0, 0, 0, 1}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({1, 2, 0, 0, 0, 4}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({1, 2, 0, 0, 0, 7}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({1, 2, 0, 0, 0, 8}));
    }

    TEST(BlockTest, FullPosition2) {
        vector<IndexList> v1{{3}, {3}};
        vector<LabelsList> v2{{WILSON_BCMUNU}, {WILSON_BCMUNU_HC}};
        auto t1 = SharedTensorData{makeVector({3}, {WILSON_BCENU}, {0.1, 0.2, 0.3}).release()};
        BlockIndexing bi{v1, v2};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        OuterElemIterator::EntryType entry0{{t1, false}};
        OuterElemIterator it0{entry0};
        auto ite0 = it0.end();
        auto chk0 = get<0>(chks[1]);
        OuterElemIterator::EntryType entry1{{t1, true}};
        OuterElemIterator it1{entry1};
        auto ite1 = it1.end();
        auto chk1 = get<0>(chks[2]);
        LabeledIndexing<AlignedIndexing> other{{3, 3}, {WILSON_BCMUNU, WILSON_BCMUNU_HC}};
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 0}));
        ++it0;
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({1, 0}));
        ++it0;
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({2, 0}));
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 0}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 1}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 2}));
    }

    TEST(BlockTest, FullPosition3) {
        vector<IndexList> v1{{2, 3}, {2, 3}};
        vector<LabelsList> v2{{FF_BD, WILSON_BCMUNU}, {FF_BD_HC, WILSON_BCMUNU_HC}};
        auto t1 = SharedTensorData{makeVector({2, 3}, {FF_BD, WILSON_BCENU}, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6}).release()};
        BlockIndexing bi{v1, v2};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        OuterElemIterator::EntryType entry0{{t1, false}};
        OuterElemIterator it0{entry0};
        auto ite0 = it0.end();
        auto chk0 = get<1>(chks[1]);
        OuterElemIterator::EntryType entry1{{t1, true}};
        OuterElemIterator it1{entry1};
        auto ite1 = it1.end();
        auto chk1 = get<1>(chks[2]);
        LabeledIndexing<AlignedIndexing> other{{2, 3, 2, 3}, {FF_BD, WILSON_BCMUNU, FF_BD_HC, WILSON_BCMUNU_HC}};
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 0, 0, 0}));
        ++it0;
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 1, 0, 0}));
        ++it0;
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 2, 0, 0}));
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 0, 0, 0}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 0, 0, 1}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 0, 0, 2}));
    }

    TEST(BlockTest, FullPosition2Bis) {
        vector<IndexList> v1{{2}, {2}};
        vector<LabelsList> v2{{WILSON_BCMUNU}, {WILSON_BCMUNU_HC}};
        auto t1 = SharedTensorData{makeVector({2}, {WILSON_BCENU}, {0.1, 0.2}).release()};
        BlockIndexing bi{v1, v2};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        OuterElemIterator::EntryType entry0{{t1, false}};
        OuterElemIterator it0{entry0};
        auto ite0 = it0.end();
        auto chk0 = get<0>(chks[1]);
        OuterElemIterator::EntryType entry1{{t1, true}};
        OuterElemIterator it1{entry1};
        auto ite1 = it1.end();
        auto chk1 = get<0>(chks[2]);
        LabeledIndexing<AlignedIndexing> other{{2, 2}, {WILSON_BCMUNU, WILSON_BCMUNU_HC}};
        EXPECT_EQ(it0.position(0), 0);
        EXPECT_FALSE(it0.isAligned(0));
        EXPECT_EQ(chk0[0], 0);
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 0}));
        ++it0;
        EXPECT_EQ(it0.position(0), 1);
        EXPECT_FALSE(it0.isAligned(0));
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({1, 0}));
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 0}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 1}));
    }

    TEST(BlockTest, FullPosition3bis) {
        vector<IndexList> v1{{3, 2}, {3, 2}};
        vector<LabelsList> v2{{FF_BD, WILSON_BCMUNU}, {FF_BD_HC, WILSON_BCMUNU_HC}};
        auto t1 = SharedTensorData{makeVector({3, 2}, {FF_BD, WILSON_BCENU}, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6}).release()};
        BlockIndexing bi{v1, v2};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        OuterElemIterator::EntryType entry0{{t1, false}};
        OuterElemIterator it0{entry0};
        auto ite0 = it0.end();
        auto chk0 = get<1>(chks[1]);
        OuterElemIterator::EntryType entry1{{t1, true}};
        OuterElemIterator it1{entry1};
        auto ite1 = it1.end();
        auto chk1 = get<1>(chks[2]);
        LabeledIndexing<AlignedIndexing> other{{3, 2, 3, 2}, {FF_BD, WILSON_BCMUNU, FF_BD_HC, WILSON_BCMUNU_HC}};
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 0, 0, 0}));
        ++it0;
        EXPECT_EQ(bi.buildFullPosition(it0, chk0), other.indicesToPos({0, 1, 0, 0}));
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 0, 0, 0}));
        ++it1;
        EXPECT_EQ(bi.buildFullPosition(it1, chk1), other.indicesToPos({0, 0, 0, 1}));
    }

    TEST(BlockTest, SplitPosition) {
        vector<IndexList> v1{{3}, {3}};
        vector<LabelsList> v2{{WILSON_BCMUNU}, {WILSON_BCMUNU_HC}};
        auto t1 = SharedTensorData{makeVector({3}, {WILSON_BCENU}, {0.1, 0.2, 0.3}).release()};
        BlockIndexing biL{v1, v2};
        vector<IndexList> v3{{2, 3}, {2, 3}};
        vector<LabelsList> v4{{FF_BD, WILSON_BCMUNU}, {FF_BD_HC, WILSON_BCMUNU_HC}};
        auto t2 = SharedTensorData{makeVector({2, 3}, {FF_BD, WILSON_BCENU}, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6}).release()};
        BlockIndexing biR{v3, v4};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        auto infoL = biL.processShifts(chks, IndexPairMember::Left);
        auto infoR = biR.processShifts(chks, IndexPairMember::Right);
        EXPECT_EQ(get<2>(infoR[0]), 2);
        OuterElemIterator::EntryType entryL0{{t1, false}};
        OuterElemIterator itL0{entryL0};
        auto itLe0 = itL0.end();
        OuterElemIterator::EntryType entryR0{{t2, false}};
        OuterElemIterator itR0{entryR0};
        auto itRe0 = itR0.end();
        IndexList inners(1);
        FlagList innerAdds(1);
        auto tmpLeft = biL.splitPosition(itL0, chks[1], get<0>(infoL[0]), get<1>(infoL[0]), inners, innerAdds);
        EXPECT_EQ(inners[0], 0);
        EXPECT_TRUE(innerAdds[0]);
        EXPECT_EQ(tmpLeft, 0);
        PositionType tmpRight;
        for (size_t i = 0; i < 6; ++i) {
            tmpRight = biR.splitPosition(itR0, chks[1], get<0>(infoR[0]), get<1>(infoR[0]), inners, innerAdds, true);
            if (i % 3 == 0) {
                EXPECT_EQ(tmpRight, i / 3);
            } else {
                EXPECT_EQ(tmpRight, numeric_limits<size_t>::max());
            }
            ++itR0;
        }
        ++itL0;
        tmpLeft = biL.splitPosition(itL0, chks[1], get<0>(infoL[0]), get<1>(infoL[0]), inners, innerAdds);
        EXPECT_EQ(inners[0], 1);
        EXPECT_TRUE(innerAdds[0]);
        EXPECT_EQ(tmpLeft, 0);
        OuterElemIterator itR0b{entryR0};
        for (size_t i = 0; i < 6; ++i) {
            tmpRight = biR.splitPosition(itR0b, chks[1], get<0>(infoR[0]), get<1>(infoR[0]), inners, innerAdds, true);
            if (i % 3 == 1) {
                EXPECT_EQ(tmpRight, i / 3);
            } else {
                EXPECT_EQ(tmpRight, numeric_limits<size_t>::max());
            }
            ++itR0b;
        }
        ++itL0;
        tmpLeft = biL.splitPosition(itL0, chks[1], get<0>(infoL[0]), get<1>(infoL[0]), inners, innerAdds);
        EXPECT_EQ(inners[0], 2);
        EXPECT_TRUE(innerAdds[0]);
        EXPECT_EQ(tmpLeft, 0);
        OuterElemIterator itR0c{entryR0};
        for (size_t i = 0; i < 6; ++i) {
            tmpRight = biR.splitPosition(itR0c, chks[1], get<0>(infoR[0]), get<1>(infoR[0]), inners, innerAdds, true);
            if (i % 3 == 2) {
                EXPECT_EQ(tmpRight, i / 3);
            } else {
                EXPECT_EQ(tmpRight, numeric_limits<size_t>::max());
            }
            ++itR0c;
        }
    }

    TEST(BlockTest, SplitPosition2) {
        vector<IndexList> v1{{3}, {3}};
        vector<LabelsList> v2{{WILSON_BCMUNU}, {WILSON_BCMUNU_HC}};
        auto t1 = SharedTensorData{makeVector({3}, {WILSON_BCENU}, {0.1, 0.2, 0.3}).release()};
        BlockIndexing biL{v1, v2};
        vector<IndexList> v3{{2, 3}, {2, 3}};
        vector<LabelsList> v4{{FF_BD, WILSON_BCMUNU}, {FF_BD_HC, WILSON_BCMUNU_HC}};
        auto t2 = SharedTensorData{makeVector({2, 3}, {FF_BD, WILSON_BCENU}, {0.1, 0.2, 0.3, 0.4, 0.5, 0.6}).release()};
        BlockIndexing biR{v3, v4};
        DotGroupList chks{DotGroupType{{}, {}, {}}, DotGroupType{{0}, {0}, {{0, 1}}}, DotGroupType{{1}, {1}, {{1, 3}}}};
        auto infoL = biL.processShifts(chks, IndexPairMember::Left);
        auto infoR = biR.processShifts(chks, IndexPairMember::Right);
        EXPECT_EQ(get<2>(infoR[1]), 2);
        OuterElemIterator::EntryType entryL1{{t1, true}};
        OuterElemIterator itL1{entryL1};
        auto itLe1 = itL1.end();
        OuterElemIterator::EntryType entryR1{{t2, true}};
        OuterElemIterator itR1{entryR1};
        auto itRe1 = itR1.end();
        IndexList inners(1);
        FlagList innerAdds(1);
        auto tmpLeft = biL.splitPosition(itL1, chks[2], get<0>(infoL[1]), get<1>(infoL[1]), inners, innerAdds);
        EXPECT_EQ(inners[0], 0);
        EXPECT_TRUE(innerAdds[0]);
        EXPECT_EQ(tmpLeft, 0);
        PositionType tmpRight;
        for (size_t i = 0; i < 6; ++i) {
            tmpRight = biR.splitPosition(itR1, chks[2], get<0>(infoR[1]), get<1>(infoR[1]), inners, innerAdds, true);
            if (i % 3 == 0) {
                EXPECT_EQ(tmpRight, i / 3);
            } else {
                EXPECT_EQ(tmpRight, numeric_limits<size_t>::max());
            }
            ++itR1;
        }
        ++itL1;
        tmpLeft = biL.splitPosition(itL1, chks[2], get<0>(infoL[1]), get<1>(infoL[1]), inners, innerAdds);
        EXPECT_EQ(inners[0], 1);
        EXPECT_TRUE(innerAdds[0]);
        EXPECT_EQ(tmpLeft, 0);
        OuterElemIterator itR1b{entryR1};
        for (size_t i = 0; i < 6; ++i) {
            tmpRight = biR.splitPosition(itR1b, chks[2], get<0>(infoR[1]), get<1>(infoR[1]), inners, innerAdds, true);
            if (i % 3 == 1) {
                EXPECT_EQ(tmpRight, i / 3);
            } else {
                EXPECT_EQ(tmpRight, numeric_limits<size_t>::max());
            }
            ++itR1b;
        }
        ++itL1;
        tmpLeft = biL.splitPosition(itL1, chks[2], get<0>(infoL[1]), get<1>(infoL[1]), inners, innerAdds);
        EXPECT_EQ(inners[0], 2);
        EXPECT_TRUE(innerAdds[0]);
        EXPECT_EQ(tmpLeft, 0);
        OuterElemIterator itR1c{entryR1};
        for (size_t i = 0; i < 6; ++i) {
            tmpRight = biR.splitPosition(itR1c, chks[2], get<0>(infoR[1]), get<1>(infoR[1]), inners, innerAdds, true);
            if (i % 3 == 2) {
                EXPECT_EQ(tmpRight, i / 3);
            } else {
                EXPECT_EQ(tmpRight, numeric_limits<size_t>::max());
            }
            ++itR1c;
        }
    }

} // namespace Hammer::MultiDimensional
