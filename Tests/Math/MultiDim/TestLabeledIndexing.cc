///
/// @file  TestLabeledIndexing.cc
/// @brief Tests for LabeledIndexing
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/SequentialIndexing.hh"
#include "Hammer/Math/MultiDim/LabeledIndexing.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    class TestableLabeledIndexing : public LabeledIndexing<SequentialIndexing> {
    public:

        TestableLabeledIndexing() = default;
        TestableLabeledIndexing(const IndexList& dimensions, const LabelsList& labels)
            : LabeledIndexing<SequentialIndexing>{dimensions, labels} {
        }
        ~TestableLabeledIndexing() = default;
        using LabeledIndexing::oppositeIndices;
        using LabeledIndexing::sameIndices;
    };


    TEST(LabeledTest, ConstructionAccess) {
        TestableLabeledIndexing t{{2, 3, 1}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX}};
        EXPECT_EQ(t.rank(), 3);
        EXPECT_EQ(t.dim(WILSON_BCENU), 2);
        EXPECT_EQ(t.dim(FF_BD), 3);
        EXPECT_EQ(t.dim(INTEGRATION_INDEX), 1);
        EXPECT_EQ(t.labels()[0], WILSON_BCENU);
        EXPECT_EQ(t.labels()[1], FF_BD);
        EXPECT_EQ(t.labels()[2], INTEGRATION_INDEX);
        TestableLabeledIndexing t2{{2, 3, 1}, {SPIN_NUTAU, SPIN_DSSD1, INTEGRATION_INDEX}};
        EXPECT_EQ(t2.dim(SPIN_NUTAU), 2);
        EXPECT_EQ(t2.dim(SPIN_DSSD1), 3);
        EXPECT_EQ(t2.dim(INTEGRATION_INDEX), 1);
        EXPECT_EQ(t2.labels()[0], SPIN_NUTAU);
        EXPECT_EQ(t2.labels()[1], SPIN_DSSD1);
        EXPECT_EQ(t2.labels()[2], INTEGRATION_INDEX);
    }

    TEST(LabeledTest, IndexRetrieval) {
        TestableLabeledIndexing t{{2, 3, 2, 3, 10},
                                  {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX}};
        auto pairs = t.oppositeIndices();
        EXPECT_EQ(pairs.size(), 2);
        EXPECT_EQ(pairs.begin()->first, WILSON_BCENU);
        EXPECT_EQ((--pairs.end())->first, FF_BD);
        EXPECT_EQ(pairs.begin()->second, WILSON_BCENU_HC);
        EXPECT_EQ((--pairs.end())->second, FF_BD_HC);
        auto indexPairs = t.getOppositeLabelPairs({{WILSON_BCENU, FF_BD}});
        EXPECT_EQ(indexPairs.size(), 2);
        EXPECT_EQ(indexPairs.begin()->first, 0);
        EXPECT_EQ((--indexPairs.end())->first, 1);
        EXPECT_EQ(indexPairs.begin()->second, 2);
        EXPECT_EQ((--indexPairs.end())->second, 3);
        TestableLabeledIndexing t2{{8, 3, 2}, {INTEGRATION_INDEX, SPIN_DSSD1, SPIN_NUTAU}};
        auto spins = t2.spinIndices();
        EXPECT_EQ(spins.size(), 2);
        EXPECT_EQ(*spins.begin(), SPIN_NUTAU);
        EXPECT_EQ(*(--spins.end()), SPIN_DSSD1);
        auto sames = t.sameIndices(t2.labels());
        EXPECT_EQ(sames.size(), 1);
        EXPECT_EQ(*sames.begin(), INTEGRATION_INDEX);
    }

    TEST(LabeledTest, IndexRetrieval2) {
        TestableLabeledIndexing t{{2, 3, 2, 3, 10},
                                  {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX}};
        TestableLabeledIndexing t2{{8, 3, 2, 2}, {INTEGRATION_INDEX, SPIN_DSSD1, SPIN_NUTAU, WILSON_BCENU}};
        EXPECT_EQ(t.labelIndex(INTEGRATION_INDEX), 4);
        EXPECT_EQ(t2.labelIndex(INTEGRATION_INDEX), 0);
        auto same = t.getSameLabelPairs(t2.labels(), {INTEGRATION_INDEX, WILSON_BCENU, FF_BD_HC, SPIN_DSSD1});
        EXPECT_EQ(same.size(), 2);
        EXPECT_EQ(same[0].first, 4);
        EXPECT_EQ(same[0].second, 0);
        EXPECT_EQ(same[1].first, 0);
        EXPECT_EQ(same[1].second, 3);
    }

    TEST(LabeledTest, SameShape) {
        TestableLabeledIndexing t{{2, 3, 2, 3, 10},
                                  {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX}};
        TestableLabeledIndexing t2{{8, 3, 2, 2}, {INTEGRATION_INDEX, SPIN_DSSD1, SPIN_NUTAU, WILSON_BCENU}};
        TestableLabeledIndexing t3{{2, 3, 2, 3, 10},
                                   {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX}};
        TestableLabeledIndexing t4{{2, 3, 2, 3, 10}, {SPIN_DSSD1, FF_BD, WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX}};
        EXPECT_TRUE(t.isSameLabelShape(t3));
        EXPECT_FALSE(t.isSameLabelShape(t2));
        EXPECT_FALSE(t.isSameLabelShape(t4));
        EXPECT_TRUE(t.isSameShape(t4));
    }

    TEST(LabeledTest, Flip) {
        TestableLabeledIndexing t{{2, 3, 2, 3}, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC}};
        TestableLabeledIndexing t2{{2, 3, 2, 3}, {WILSON_BCENU_HC, FF_BD_HC, WILSON_BCENU, FF_BD}};
        EXPECT_FALSE(t.isSameLabelShape(t2));
        t.flipLabels();
        EXPECT_TRUE(t.isSameLabelShape(t2));
    }

    TEST(LabeledTest, AddAt) {
        TestableLabeledIndexing t{{2, 3, 2, 3}, {WILSON_BCENU_HC, FF_BD_HC, WILSON_BCENU, FF_BD}};
        TestableLabeledIndexing t2{{2, 3, 2, 3, 11},
                                   {WILSON_BCENU_HC, FF_BD_HC, WILSON_BCENU, FF_BD, INTEGRATION_INDEX}};
        TestableLabeledIndexing t3{{4, 3, 4, 3, 11},
                                   {WILSON_BCENU_HC, FF_BD_HC, WILSON_BCENU, FF_BD, INTEGRATION_INDEX}};
        TestableLabeledIndexing t4{{2, 3, 11, 2, 3},
                                   {WILSON_BCENU_HC, FF_BD_HC, INTEGRATION_INDEX, WILSON_BCENU, FF_BD}};
        TestableLabeledIndexing t5{{2, 3, 11, 2, 3},
                                   {WILSON_BCMUNU_HC, FF_BD_HC, INTEGRATION_INDEX, WILSON_BCMUNU, FF_BD}};
        EXPECT_TRUE(t2.canAddAt(t.labels(), t.dims(), INTEGRATION_INDEX, 7));
        EXPECT_FALSE(t2.canAddAt(t.labels(), t.dims(), INTEGRATION_INDEX, 11));
        EXPECT_FALSE(t3.canAddAt(t.labels(), t.dims(), INTEGRATION_INDEX, 7));
        EXPECT_TRUE(t4.canAddAt(t.labels(), t.dims(), INTEGRATION_INDEX, 7));
        EXPECT_FALSE(t5.canAddAt(t.labels(), t.dims(), INTEGRATION_INDEX, 7));
    }

} // namespace Hammer::MultiDimensional
