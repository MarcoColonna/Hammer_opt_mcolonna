///
/// @file  TestAlignedIndexing.cc
/// @brief Tests for AlignedIndexing
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define ALIGNEDINDEXING_FRIENDS                   \
    FRIEND_TEST(AlignedTest, ConstructionAccess); \
    FRIEND_TEST(AlignedTest, MaxIndex);           \
    FRIEND_TEST(AlignedTest, Pads);               \
    FRIEND_TEST(AlignedTest, Pads2);              \
    FRIEND_TEST(AlignedTest, Pads3);              \
    FRIEND_TEST(AlignedTest, ValidIndices);       \
    FRIEND_TEST(AlignedTest, PositionConversion); \
    FRIEND_TEST(AlignedTest, Shifts);             \
    FRIEND_TEST(AlignedTest, Splits);             \
    FRIEND_TEST(AlignedTest, Extend)


#include "Hammer/Math/MultiDim/AlignedIndexing.hh"

using namespace std;

namespace Hammer::MultiDimensional {


    TEST(AlignedTest, ConstructionAccess) {
        AlignedIndexing t{{2, 3, 1}};
        EXPECT_EQ(t.rank(), 3);
        EXPECT_EQ(t.dim(0), 2);
        EXPECT_EQ(t.dim(1), 3);
        EXPECT_EQ(t.dim(2), 1);
        EXPECT_EQ(t.dims()[0], 2);
        EXPECT_EQ(t.dims()[1], 3);
        EXPECT_EQ(t.dims()[2], 1);
        EXPECT_EQ(t.numValues(), 6);
    }

    TEST(AlignedTest, MaxIndex) {
        AlignedIndexing t{{2, 3, 1}};
        EXPECT_EQ(t.maxIndex(false), 2 * 3 - 1);
        EXPECT_EQ(t.maxIndex(), 2 * 4 - 1);
    }

    TEST(AlignedTest, Pads) {
        AlignedIndexing t{{2, 3, 1}};
        EXPECT_EQ(t.indicesToPos({0, 0, 0}), 0);
        EXPECT_EQ(t.indicesToPos({0, 1, 0}), 1);
        EXPECT_EQ(t.indicesToPos({0, 2, 0}), 2);
        EXPECT_EQ(t.indicesToPos({1, 0, 0}), (1 << 2) + 0);
        EXPECT_EQ(t.indicesToPos({1, 1, 0}), (1 << 2) + 1);
        EXPECT_EQ(t.indicesToPos({1, 2, 0}), (1 << 2) + 2);
        IndexList result(3);
        t.posToIndices(0, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(1, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(2, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\2');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices((1 << 2) + 0, result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices((1 << 2) + 1, result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices((1 << 2) + 2, result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\2');
        EXPECT_EQ(result[2], '\0');
        EXPECT_EQ(t.ithIndexInPos((1 << 2) + 2, '\0'), '\1');
        EXPECT_EQ(t.ithIndexInPos((1 << 2) + 2, '\1'), '\2');
        EXPECT_EQ(t.ithIndexInPos((1 << 2) + 2, '\2'), '\0');
    }

    TEST(AlignedTest, Pads2) {
        AlignedIndexing t{{7, 3, 4}};
        EXPECT_EQ(t.indicesToPos({0, 0, 0}), 0);
        EXPECT_EQ(t.indicesToPos({0, 0, 3}), 3);
        EXPECT_EQ(t.indicesToPos({0, 1, 0}), 1 << 2);
        EXPECT_EQ(t.indicesToPos({0, 2, 0}), 2 << 2);
        EXPECT_EQ(t.indicesToPos({1, 0, 0}), 1 << (2 + 2));
        EXPECT_EQ(t.indicesToPos({1, 1, 0}), (1 << (2 + 2)) + (1 << 2));
        EXPECT_EQ(t.indicesToPos({1, 2, 0}), (1 << (2 + 2)) + (2 << 2));
        IndexList result(3);
        t.posToIndices(0, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(1 << 2, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(2 << 2, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\2');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(1 << (2 + 2), result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices((1 << (2 + 2)) + (1 << 2), result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices((1 << (2 + 2)) + (2 << 2), result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\2');
        EXPECT_EQ(result[2], '\0');
        EXPECT_EQ(t.ithIndexInPos((1 << (2 + 2)) + (2 << 2), '\0'), '\1');
        EXPECT_EQ(t.ithIndexInPos((1 << (2 + 2)) + (2 << 2), '\1'), '\2');
        EXPECT_EQ(t.ithIndexInPos((1 << (2 + 2)) + (2 << 2), '\2'), '\0');
    }

    TEST(AlignedTest, Pads3) {
        AlignedIndexing t{{7, 3, 2}};
        EXPECT_EQ(t.indicesToPos({0, 0, 0}), 0);
        EXPECT_EQ(t.indicesToPos({0, 0, 1}), 1);
        EXPECT_EQ(t.indicesToPos({0, 1, 0}), 1 << 1);
        EXPECT_EQ(t.indicesToPos({0, 2, 0}), 2 << 1);
        EXPECT_EQ(t.indicesToPos({1, 0, 0}), 1 << (2 + 1));
        EXPECT_EQ(t.indicesToPos({1, 1, 0}), (1 << (2 + 1)) + (1 << 1));
        EXPECT_EQ(t.indicesToPos({1, 2, 0}), (1 << (2 + 1)) + (2 << 1));
    }

    TEST(AlignedTest, ValidIndices) {
        AlignedIndexing t{{2, 3, 1}};
        EXPECT_TRUE(t.checkValidIndices({1, 2, 0}));
        EXPECT_FALSE(t.checkValidIndices({2, 2, 0}));
        EXPECT_FALSE(t.checkValidIndices({1, 3, 0}));
        EXPECT_FALSE(t.checkValidIndices({1, 2, 1}));
        EXPECT_FALSE(t.checkValidIndices({2, 2, 1}));
    }

    TEST(AlignedTest, PositionConversion) {
        AlignedIndexing t{{2, 7, 3}};
        EXPECT_EQ(t.alignedPosToPos(t.indicesToPos({1, 6, 2})), 2 + 6 * 3 + 1 * 7 * 3);
        EXPECT_EQ(t.alignedPosToPos(t.indicesToPos({0, 2, 1})), 1 + 2 * 3 + 0 * 7 * 3);
        EXPECT_EQ(t.posToAlignedPos(2 + 6 * 3 + 1 * 7 * 3), t.indicesToPos({1, 6, 2}));
        EXPECT_EQ(t.posToAlignedPos(1 + 2 * 3 + 0 * 7 * 3), t.indicesToPos({0, 2, 1}));
    }

    TEST(AlignedTest, Shifts) {
        AlignedIndexing t1{{2, 7, 3}};
        AlignedIndexing t2{{11, 2, 3, 4}};
        AlignedIndexing t3{{11, 2, 11, 3, 2}};
        auto res1 = t1.processShifts({{0, 1}}, IndexPairMember::Left);
        auto res2 = t2.processShifts({{0, 1}}, IndexPairMember::Right);
        auto res3 = t3.processShifts({{0, 2}}, IndexPairMember::Both);
        EXPECT_EQ(get<0>(res1).size(), 3);
        EXPECT_EQ(get<1>(res1).size(), 3);
        EXPECT_EQ(get<0>(res2).size(), 4);
        EXPECT_EQ(get<1>(res2).size(), 4);
        EXPECT_EQ(get<0>(res3).size(), 5);
        EXPECT_EQ(get<1>(res3).size(), 5);
        // t1
        // contracted
        EXPECT_EQ(get<0>(res1)[0], 0);
        EXPECT_FALSE(get<1>(res1)[0]);
        // uncontracted
        EXPECT_EQ(get<0>(res1)[1], 2);
        EXPECT_TRUE(get<1>(res1)[1]);
        // uncontracted
        EXPECT_EQ(get<0>(res1)[2], 0);
        EXPECT_TRUE(get<1>(res1)[2]);
        // t2
        // uncontracted
        EXPECT_EQ(get<0>(res2)[0], 4);
        EXPECT_TRUE(get<1>(res2)[0]);
        // contracted
        EXPECT_EQ(get<0>(res2)[1], 0);
        EXPECT_FALSE(get<1>(res2)[1]);
        // uncontracted
        EXPECT_EQ(get<0>(res2)[2], 2);
        EXPECT_TRUE(get<1>(res2)[2]);
        // uncontracted
        EXPECT_EQ(get<0>(res2)[3], 0);
        EXPECT_TRUE(get<1>(res2)[3]);
        // t3
        // contracted
        EXPECT_EQ(get<0>(res3)[0], 0);
        EXPECT_FALSE(get<1>(res3)[0]);
        // uncontracted
        EXPECT_EQ(get<0>(res3)[1], 3);
        EXPECT_TRUE(get<1>(res3)[1]);
        // contracted
        EXPECT_EQ(get<0>(res3)[2], 0);
        EXPECT_FALSE(get<1>(res3)[2]);
        // uncontracted
        EXPECT_EQ(get<0>(res3)[3], 1);
        EXPECT_TRUE(get<1>(res3)[3]);
        // uncontracted
        EXPECT_EQ(get<0>(res3)[4], 0);
        EXPECT_TRUE(get<1>(res3)[4]);
        res1 = t1.processShifts({{2, 2}, {0, 1}}, IndexPairMember::Left);
        res2 = t2.processShifts({{0, 1}, {2, 2}}, IndexPairMember::Right);
        res3 = t3.processShifts({{0, 2}, {1, 4}}, IndexPairMember::Both);
        // t1
        // contracted
        EXPECT_EQ(get<0>(res1)[0], 1);
        EXPECT_FALSE(get<1>(res1)[0]);
        // uncontracted
        EXPECT_EQ(get<0>(res1)[1], 0);
        EXPECT_TRUE(get<1>(res1)[1]);
        // contracted
        EXPECT_EQ(get<0>(res1)[2], 0);
        EXPECT_FALSE(get<1>(res1)[2]);
        // t2
        // uncontracted
        EXPECT_EQ(get<0>(res2)[0], 2);
        EXPECT_TRUE(get<1>(res2)[0]);
        // contracted
        EXPECT_EQ(get<0>(res2)[1], 0);
        EXPECT_FALSE(get<1>(res2)[1]);
        // contracted
        EXPECT_EQ(get<0>(res2)[2], 1);
        EXPECT_FALSE(get<1>(res2)[2]);
        // uncontracted
        EXPECT_EQ(get<0>(res2)[3], 0);
        EXPECT_TRUE(get<1>(res2)[3]);
        // t3
        // contracted
        EXPECT_EQ(get<0>(res3)[0], 0);
        EXPECT_FALSE(get<1>(res3)[0]);
        // contracted
        EXPECT_EQ(get<0>(res3)[1], 1);
        EXPECT_FALSE(get<1>(res3)[1]);
        // contracted
        EXPECT_EQ(get<0>(res3)[2], 0);
        EXPECT_FALSE(get<1>(res3)[2]);
        // uncontracted
        EXPECT_EQ(get<0>(res3)[3], 0);
        EXPECT_TRUE(get<1>(res3)[3]);
        // contracted
        EXPECT_EQ(get<0>(res3)[4], 1);
        EXPECT_FALSE(get<1>(res3)[4]);
    }

    TEST(AlignedTest, Splits) {
        AlignedIndexing t1{{2, 7, 3}};
        AlignedIndexing t2{{11, 2, 3, 4}};
        AlignedIndexing t3{{11, 2, 11, 3, 2}};
        auto res1 = t1.processShifts({{0, 1}}, IndexPairMember::Left);
        auto res2 = t2.processShifts({{0, 1}}, IndexPairMember::Right);
        auto res3 = t3.processShifts({{0, 2}}, IndexPairMember::Both);
        IndexList innerA = {0};
        IndexList innerB = {0, 0};
        IndexList innerC = {0, 0, 0, 0};
        vector<bool> innerAddA = {false};
        vector<bool> innerAddB = {false, false};
        vector<bool> innerAddC = {false, false, false, false};
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto split1 =
                        t1.splitPosition(t1.indicesToPos({i1, i2, i3}), get<0>(res1), get<1>(res1), innerA, innerAddA);
                    EXPECT_EQ(split1, (i2 << 2) + i3);
                    EXPECT_EQ(innerA.size(), 1);
                    EXPECT_EQ(innerA[0], i1);
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 4; ++i4) {
                        auto split2 = t2.splitPosition(t2.indicesToPos({i1, i2, i3, i4}), get<0>(res2), get<1>(res2),
                                                       innerA, innerAddA);
                        EXPECT_EQ(split2, (i1 << 4) + (i3 << 2) + i4);
                        EXPECT_EQ(innerA.size(), 1);
                        EXPECT_EQ(innerA[0], i2);
                    }
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 11; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        for (IndexType i5 = 0; i5 < 2; ++i5) {
                            auto split3 = t3.splitPosition(t3.indicesToPos({i1, i2, i3, i4, i5}), get<0>(res3),
                                                           get<1>(res3), innerB, innerAddB);
                            if (i1 == i3) {
                                EXPECT_EQ(split3, (i2 << 3) + (i4 << 1) + i5);
                                EXPECT_EQ(innerB.size(), 2);
                                EXPECT_EQ(innerB[0], i1);
                                EXPECT_EQ(innerB[1], 0);
                            } else {
                                EXPECT_EQ(split3, numeric_limits<size_t>::max());
                            }
                        }
                    }
                }
            }
        }
        res1 = t1.processShifts({{2, 2}, {0, 1}}, IndexPairMember::Left);
        res2 = t2.processShifts({{2, 2}, {0, 1}}, IndexPairMember::Right);
        res3 = t3.processShifts({{0, 2}, {1, 4}}, IndexPairMember::Both);
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto split1 =
                        t1.splitPosition(t1.indicesToPos({i1, i2, i3}), get<0>(res1), get<1>(res1), innerB, innerAddB);
                    EXPECT_EQ(split1, i2);
                    EXPECT_EQ(innerB.size(), 2);
                    EXPECT_EQ(innerB[0], i3);
                    EXPECT_EQ(innerB[1], i1);
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 4; ++i4) {
                        auto split2 = t2.splitPosition(t2.indicesToPos({i1, i2, i3, i4}), get<0>(res2), get<1>(res2),
                                                       innerB, innerAddB);
                        EXPECT_EQ(split2, (i1 << 2) + i4);
                        EXPECT_EQ(innerB.size(), 2);
                        EXPECT_EQ(innerB[0], i3);
                        EXPECT_EQ(innerB[1], i2);
                    }
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 11; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        for (IndexType i5 = 0; i5 < 2; ++i5) {
                            auto split3 = t3.splitPosition(t3.indicesToPos({i1, i2, i3, i4, i5}), get<0>(res3),
                                                           get<1>(res3), innerC, innerAddC);
                            if (i1 == i3 && i2 == i5) {
                                EXPECT_EQ(split3, i4);
                                EXPECT_EQ(innerC.size(), 4);
                                EXPECT_EQ(innerC[0], i1);
                                EXPECT_EQ(innerC[1], i2);
                                EXPECT_EQ(innerC[2], 0);
                                EXPECT_EQ(innerC[3], 0);
                            } else {
                                EXPECT_EQ(split3, numeric_limits<size_t>::max());
                            }
                        }
                    }
                }
            }
        }
        IndexList vec1 = {1, 1};
        IndexList vec2 = {2, 0};
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 4; ++i4) {
                        auto split2a = t2.splitPosition(t2.indicesToPos({i1, i2, i3, i4}), get<0>(res2), get<1>(res2),
                                                        vec1, innerAddB, true);
                        auto split2b = t2.splitPosition(t2.indicesToPos({i1, i2, i3, i4}), get<0>(res2), get<1>(res2),
                                                        vec2, innerAddB, true);
                        if (i3 != 1 || i2 != 1) {
                            EXPECT_EQ(split2a, numeric_limits<size_t>::max());
                        } else {
                            EXPECT_EQ(split2a, (i1 << 2) + i4);
                        }
                        if (i3 != 2 || i2 != 0) {
                            EXPECT_EQ(split2b, numeric_limits<size_t>::max());
                        } else {
                            EXPECT_EQ(split2b, (i1 << 2) + i4);
                        }
                    }
                }
            }
        }
    }

    TEST(AlignedTest, Extend) {
        AlignedIndexing t1{{2, 7, 3}};
        AlignedIndexing t2{{2, 7, 11, 3}};
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 11; ++i4) {
                        EXPECT_EQ(t2.extendAlignedPosition(t1.indicesToPos({i1, i2, i3}), 2, i4),
                                  t2.indicesToPos({i1, i2, i4, i3}));
                        EXPECT_EQ(t2.extendPosition(t1.alignedPosToPos(t1.indicesToPos({i1, i2, i3})), 3, 2, i4),
                                  t2.indicesToPos({i1, i2, i4, i3}));
                    }
                }
            }
        }
    }

} // namespace Hammer::MultiDimensional
