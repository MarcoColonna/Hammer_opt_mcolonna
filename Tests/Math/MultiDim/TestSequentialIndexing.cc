///
/// @file  TestSequentialIndexing.cc
/// @brief Tests for SequentialIndexing
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define SEQUENTIALINDEXING_FRIENDS                   \
    FRIEND_TEST(SequentialTest, ConstructionAccess); \
    FRIEND_TEST(SequentialTest, Strides);            \
    FRIEND_TEST(SequentialTest, Strides2);           \
    FRIEND_TEST(SequentialTest, Pads3);              \
    FRIEND_TEST(SequentialTest, ValidIndices);       \
    FRIEND_TEST(SequentialTest, Shifts);             \
    FRIEND_TEST(SequentialTest, Splits);             \
    FRIEND_TEST(SequentialTest, Extend)

#include "Hammer/Math/MultiDim/SequentialIndexing.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    TEST(SequentialTest, ConstructionAccess) {
        SequentialIndexing t{{2, 3, 1}};
        EXPECT_EQ(t.rank(), 3);
        EXPECT_EQ(t.dim(0), 2);
        EXPECT_EQ(t.dim(1), 3);
        EXPECT_EQ(t.dim(2), 1);
        EXPECT_EQ(t.dims()[0], 2);
        EXPECT_EQ(t.dims()[1], 3);
        EXPECT_EQ(t.dims()[2], 1);
        EXPECT_EQ(t.numValues(), 6);
    }

    TEST(SequentialTest, Strides) {
        SequentialIndexing t{{2, 3, 1}};
        EXPECT_EQ(t.indicesToPos({0, 0, 0}), 0);
        EXPECT_EQ(t.indicesToPos({0, 1, 0}), 1);
        EXPECT_EQ(t.indicesToPos({0, 2, 0}), 2 * 1);
        EXPECT_EQ(t.indicesToPos({1, 0, 0}), 3);
        EXPECT_EQ(t.indicesToPos({1, 1, 0}), 3 + 1);
        EXPECT_EQ(t.indicesToPos({1, 2, 0}), 3 + 2 * 1);
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
        t.posToIndices(3, result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(4, result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(5, result);
        EXPECT_EQ(result[0], '\1');
        EXPECT_EQ(result[1], '\2');
        EXPECT_EQ(result[2], '\0');
        EXPECT_EQ(t.ithIndexInPos(5, '\0'), '\1');
        EXPECT_EQ(t.ithIndexInPos(5, '\1'), '\2');
        EXPECT_EQ(t.ithIndexInPos(5, '\2'), '\0');
    }


    TEST(SequentialTest, Strides2) {
        SequentialIndexing t{{7, 3, 4}};
        EXPECT_EQ(t.indicesToPos({0, 0, 0}), 0);
        EXPECT_EQ(t.indicesToPos({0, 0, 3}), 3);
        EXPECT_EQ(t.indicesToPos({0, 1, 0}), 1 * 4 + 0);
        EXPECT_EQ(t.indicesToPos({0, 2, 0}), 2 * 4 + 0);
        EXPECT_EQ(t.indicesToPos({1, 0, 0}), 1 * (3 * 4) + 0 * 4 + 0);
        EXPECT_EQ(t.indicesToPos({1, 1, 0}), 1 * (3 * 4) + 1 * 4 + 0);
        EXPECT_EQ(t.indicesToPos({1, 2, 0}), 1 * (3 * 4) + 2 * 4 + 0);
        IndexList result(3);
        t.posToIndices(0, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(1, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\1');
        t.posToIndices(2, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\2');
        t.posToIndices(3, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\3');
        t.posToIndices(4, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\0');
        t.posToIndices(5, result);
        EXPECT_EQ(result[0], '\0');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\1');
        t.posToIndices(25, result);
        EXPECT_EQ(result[0], '\2');
        EXPECT_EQ(result[1], '\0');
        EXPECT_EQ(result[2], '\1');
        t.posToIndices(79, result);
        EXPECT_EQ(result[0], '\6');
        EXPECT_EQ(result[1], '\1');
        EXPECT_EQ(result[2], '\3');
        EXPECT_EQ(t.ithIndexInPos(5, '\0'), '\0');
        EXPECT_EQ(t.ithIndexInPos(5, '\1'), '\1');
        EXPECT_EQ(t.ithIndexInPos(5, '\2'), '\1');
        EXPECT_EQ(t.ithIndexInPos(79, '\0'), '\6');
        EXPECT_EQ(t.ithIndexInPos(79, '\1'), '\1');
        EXPECT_EQ(t.ithIndexInPos(79, '\2'), '\3');
    }

    TEST(SequentialTest, Pads3) {
        SequentialIndexing t{{7, 3, 2}};
        EXPECT_EQ(t.indicesToPos({0, 0, 0}), 0);
        EXPECT_EQ(t.indicesToPos({0, 0, 1}), 1);
        EXPECT_EQ(t.indicesToPos({0, 1, 0}), 1 * 2);
        EXPECT_EQ(t.indicesToPos({0, 2, 0}), 2 * 2);
        EXPECT_EQ(t.indicesToPos({1, 0, 0}), 1 * 6);
        EXPECT_EQ(t.indicesToPos({1, 1, 0}), (1 * 6) + (1 * 2));
        EXPECT_EQ(t.indicesToPos({1, 2, 0}), (1 * 6) + (2 * 2));
        EXPECT_EQ(t.stride(0), 6);
        EXPECT_EQ(t.stride(1), 2);
        EXPECT_EQ(t.stride(2), 1);
    }

    TEST(SequentialTest, ValidIndices) {
        SequentialIndexing t{{2, 3, 1}};
        EXPECT_TRUE(t.checkValidIndices({1, 2, 0}));
        EXPECT_FALSE(t.checkValidIndices({2, 2, 0}));
        EXPECT_FALSE(t.checkValidIndices({1, 3, 0}));
        EXPECT_FALSE(t.checkValidIndices({1, 2, 1}));
        EXPECT_FALSE(t.checkValidIndices({2, 2, 1}));
    }

    TEST(SequentialTest, Shifts) {
        IndexPairList v1 = {{0, 1}};
        IndexPairList v2 = {{0, 1}, {2, 2}};
        IndexPairList v3 = {{0, 2}};
        IndexPairList v4 = {{0, 2}, {1, 4}};
        SequentialIndexing t1{{2, 7, 3}};
        SequentialIndexing t2{{11, 2, 3, 4}};
        SequentialIndexing t3{{11, 2, 11, 3, 2}};
        auto res1a = t1.getInnerOuterStrides(v1, t1.strides());
        auto res1b = t1.getInnerOuterStrides(v1, t2.strides());
        auto res2a = t2.getInnerOuterStrides(v1, t2.strides());
        auto res2b = t2.getOuterStrides2nd(v1);
        auto res3 = t3.getInnerOuterStrides(v3, t3.strides(), true);
        EXPECT_EQ(res1a.size(), 3);
        EXPECT_EQ(res1b.size(), 3);
        EXPECT_EQ(res2a.size(), 4);
        EXPECT_EQ(res2b.size(), 4 - 1);
        EXPECT_EQ(res3.size(), 5);
        // res1a
        EXPECT_EQ(get<0>(res1a[0]), 21);
        EXPECT_EQ(get<1>(res1a[0]), 3);
        EXPECT_TRUE(get<2>(res1a[0]));
        EXPECT_EQ(get<0>(res1a[1]), 0);
        EXPECT_EQ(get<1>(res1a[1]), 0);
        EXPECT_FALSE(get<2>(res1a[1]));
        EXPECT_EQ(get<0>(res1a[2]), 1);
        EXPECT_EQ(get<1>(res1a[2]), 1);
        EXPECT_FALSE(get<2>(res1a[2]));
        // res1b
        EXPECT_EQ(get<0>(res1b[0]), 21);
        EXPECT_EQ(get<1>(res1b[0]), 12);
        EXPECT_TRUE(get<2>(res1b[0]));
        EXPECT_EQ(get<0>(res1b[1]), 0);
        EXPECT_EQ(get<1>(res1b[1]), 0);
        EXPECT_FALSE(get<2>(res1b[1]));
        EXPECT_EQ(get<0>(res1b[2]), 1);
        EXPECT_EQ(get<1>(res1b[2]), 1);
        EXPECT_FALSE(get<2>(res1b[2]));
        // res2a
        EXPECT_EQ(get<0>(res2a[0]), 24);
        EXPECT_EQ(get<1>(res2a[0]), 12);
        EXPECT_TRUE(get<2>(res2a[0]));
        EXPECT_EQ(get<0>(res2a[1]), 0);
        EXPECT_EQ(get<1>(res2a[1]), 0);
        EXPECT_FALSE(get<2>(res2a[1]));
        EXPECT_EQ(get<0>(res2a[2]), 0);
        EXPECT_EQ(get<1>(res2a[2]), 0);
        EXPECT_FALSE(get<2>(res2a[2]));
        EXPECT_EQ(get<0>(res2a[3]), 1);
        EXPECT_EQ(get<1>(res2a[3]), 1);
        EXPECT_FALSE(get<2>(res2a[3]));
        // res2b
        EXPECT_EQ(res2b[0].first, 3 * 4);
        EXPECT_EQ(res2b[0].second, 2 * 3 * 4);
        EXPECT_EQ(res2b[1].first, 4);
        EXPECT_EQ(res2b[1].second, 4);
        EXPECT_EQ(res2b[2].first, 1);
        EXPECT_EQ(res2b[2].second, 1);
        // res3
        EXPECT_EQ(get<0>(res3[0]), 132);
        EXPECT_EQ(get<1>(res3[0]), 6);
        EXPECT_TRUE(get<2>(res3[0]));
        EXPECT_EQ(get<0>(res3[1]), 66);
        EXPECT_EQ(get<1>(res3[1]), 6);
        EXPECT_FALSE(get<2>(res3[1]));
        EXPECT_EQ(get<0>(res3[2]), 6);
        EXPECT_EQ(get<1>(res3[2]), -6);
        EXPECT_TRUE(get<2>(res3[2]));
        EXPECT_EQ(get<0>(res3[3]), 0);
        EXPECT_EQ(get<1>(res3[3]), 0);
        EXPECT_FALSE(get<2>(res3[3]));
        EXPECT_EQ(get<0>(res3[4]), 1);
        EXPECT_EQ(get<1>(res3[4]), 1);
        EXPECT_FALSE(get<2>(res3[4]));
        res1a = t1.getInnerOuterStrides(v2, t1.strides());
        res1b = t1.getInnerOuterStrides(v2, t2.strides());
        res2a = t2.getInnerOuterStrides(v2, t2.strides());
        res2b = t2.getOuterStrides2nd(v2);
        res3 = t3.getInnerOuterStrides(v4, t3.strides(), true);
        EXPECT_EQ(res1a.size(), 3);
        EXPECT_EQ(res1b.size(), 3);
        EXPECT_EQ(res2a.size(), 4);
        EXPECT_EQ(res2b.size(), 4 - 2);
        EXPECT_EQ(res3.size(), 5);
        // res1a
        EXPECT_EQ(get<0>(res1a[0]), 21);
        EXPECT_EQ(get<1>(res1a[0]), 3);
        EXPECT_TRUE(get<2>(res1a[0]));
        EXPECT_EQ(get<0>(res1a[1]), 3);
        EXPECT_EQ(get<1>(res1a[1]), 1);
        EXPECT_FALSE(get<2>(res1a[1]));
        EXPECT_EQ(get<0>(res1a[2]), 1);
        EXPECT_EQ(get<1>(res1a[2]), 1);
        EXPECT_TRUE(get<2>(res1a[2]));
        // res1b
        EXPECT_EQ(get<0>(res1b[0]), 21);
        EXPECT_EQ(get<1>(res1b[0]), 12);
        EXPECT_TRUE(get<2>(res1b[0]));
        EXPECT_EQ(get<0>(res1b[1]), 3);
        EXPECT_EQ(get<1>(res1b[1]), 1);
        EXPECT_FALSE(get<2>(res1b[1]));
        EXPECT_EQ(get<0>(res1b[2]), 1);
        EXPECT_EQ(get<1>(res1b[2]), 4);
        EXPECT_TRUE(get<2>(res1b[2]));
        // res2a
        EXPECT_EQ(get<0>(res2a[0]), 24);
        EXPECT_EQ(get<1>(res2a[0]), 12);
        EXPECT_TRUE(get<2>(res2a[0]));
        EXPECT_EQ(get<0>(res2a[1]), 12);
        EXPECT_EQ(get<1>(res2a[1]), 4);
        EXPECT_FALSE(get<2>(res2a[1]));
        EXPECT_EQ(get<0>(res2a[2]), 4);
        EXPECT_EQ(get<1>(res2a[2]), 4);
        EXPECT_TRUE(get<2>(res2a[2]));
        EXPECT_EQ(get<0>(res2a[3]), 1);
        EXPECT_EQ(get<1>(res2a[3]), 1);
        EXPECT_FALSE(get<2>(res2a[3]));
        // // res2b
        EXPECT_EQ(res2b[0].first, 4);
        EXPECT_EQ(res2b[0].second, 2 * 3 * 4);
        EXPECT_EQ(res2b[1].first, 1);
        EXPECT_EQ(res2b[1].second, 1);
        // // res3
        EXPECT_EQ(get<0>(res3[0]), 132);
        EXPECT_EQ(get<1>(res3[0]), 6);
        EXPECT_TRUE(get<2>(res3[0]));
        EXPECT_EQ(get<0>(res3[1]), 66);
        EXPECT_EQ(get<1>(res3[1]), 1);
        EXPECT_TRUE(get<2>(res3[1]));
        EXPECT_EQ(get<0>(res3[2]), 6);
        EXPECT_EQ(get<1>(res3[2]), -6);
        EXPECT_TRUE(get<2>(res3[2]));
        EXPECT_EQ(get<0>(res3[3]), 2);
        EXPECT_EQ(get<1>(res3[3]), 1);
        EXPECT_FALSE(get<2>(res3[3]));
        EXPECT_EQ(get<0>(res3[4]), 1);
        EXPECT_EQ(get<1>(res3[4]), -1);
        EXPECT_TRUE(get<2>(res3[4]));
    }

    TEST(SequentialTest, Splits) {
        IndexPairList v1 = {{0, 1}};
        IndexPairList v2 = {{0, 1}, {2, 2}};
        IndexPairList v3 = {{0, 2}};
        IndexPairList v4 = {{0, 2}, {1, 4}};
        SequentialIndexing t1{{2, 7, 3}};
        SequentialIndexing t2{{11, 2, 3, 4}};
        SequentialIndexing t3{{11, 2, 11, 3, 2}};
        auto res1 = t1.getInnerOuterStrides(v1, t2.strides());
        auto res2b = t2.getOuterStrides2nd(v1);
        auto res3 = t3.getInnerOuterStrides(v3, t3.strides(), true);
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto split1 = SequentialIndexing::splitPosition(t1.indicesToPos({i1, i2, i3}), res1);
                    EXPECT_EQ(split1.first, i2 * 3 + i3);
                    EXPECT_EQ(split1.second, i1 * 12);
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 4; ++i4) {
                        EXPECT_EQ(
                            t2.build2ndPosition(i1 * 12u + i3 * 4u + i4, static_cast<PositionType>(0 * 12), res2b),
                            i1 * 24 + 0 * 12 + i3 * 4 + i4);
                        EXPECT_EQ(
                            t2.build2ndPosition(i1 * 12u + i3 * 4u + i4, static_cast<PositionType>(1 * 12), res2b),
                            i1 * 24 + 1u * 12 + i3 * 4 + i4);
                    }
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 11; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        for (IndexType i5 = 0; i5 < 2; ++i5) {
                            auto split3 =
                                SequentialIndexing::splitPosition(t3.indicesToPos({i1, i2, i3, i4, i5}), res3);
                            EXPECT_EQ(split3.first, i2 * 6 + i4 * 2 + i5);
                            if (i1 == i3) {
                                EXPECT_EQ(split3.second, 0);
                            } else {
                                EXPECT_NE(split3.second, 0);
                            }
                        }
                    }
                }
            }
        }
        EXPECT_EQ(res1.size(), 3);
        EXPECT_EQ(res3.size(), 5);
        res1 = t1.getInnerOuterStrides(v2, t2.strides());
        res2b = t2.getOuterStrides2nd(v2);
        res3 = t3.getInnerOuterStrides(v4, t3.strides(), true);
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto split1 = SequentialIndexing::splitPosition(t1.indicesToPos({i1, i2, i3}), res1);
                    EXPECT_EQ(split1.first, i2);
                    EXPECT_EQ(split1.second, i1 * 12 + i3 * 4);
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 4; ++i4) {
                        EXPECT_EQ(t2.build2ndPosition(i1 * 4u + i4, 0 * 12 + 0 * 4, res2b),
                                  i1 * 24 + 0 * 12 + 0 * 4 + i4);
                        EXPECT_EQ(t2.build2ndPosition(i1 * 4u + i4, 1 * 12 + 0 * 4, res2b),
                                  i1 * 24 + 1 * 12 + 0 * 4 + i4);
                        EXPECT_EQ(t2.build2ndPosition(i1 * 4u + i4, 0 * 12 + 1 * 4, res2b),
                                  i1 * 24 + 0 * 12 + 1 * 4 + i4);
                        EXPECT_EQ(t2.build2ndPosition(i1 * 4u + i4, 1 * 12 + 1 * 4, res2b),
                                  i1 * 24 + 1 * 12 + 1 * 4 + i4);
                        EXPECT_EQ(t2.build2ndPosition(i1 * 4u + i4, 0 * 12 + 2 * 4, res2b),
                                  i1 * 24 + 0 * 12 + 2 * 4 + i4);
                        EXPECT_EQ(t2.build2ndPosition(i1 * 4u + i4, 1 * 12 + 2 * 4, res2b),
                                  i1 * 24 + 1 * 12 + 2 * 4 + i4);
                    }
                }
            }
        }
        for (IndexType i1 = 0; i1 < 11; ++i1) {
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                for (IndexType i3 = 0; i3 < 11; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        for (IndexType i5 = 0; i5 < 2; ++i5) {
                            auto split3 =
                                SequentialIndexing::splitPosition(t3.indicesToPos({i1, i2, i3, i4, i5}), res3);
                            EXPECT_EQ(split3.first, i4);
                            if (i1 == i3 && i2 == i5) {
                                EXPECT_EQ(split3.second, 0);
                            } else {
                                EXPECT_NE(split3.second, 0);
                            }
                        }
                    }
                }
            }
        }
    }

    TEST(SequentialTest, Extend) {
        SequentialIndexing t1{{2, 7, 3}};
        SequentialIndexing t2{{2, 7, 11, 3}};
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 7; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 11; ++i4) {
                        EXPECT_EQ(t2.extendPosition(t1.indicesToPos({i1, i2, i3}), 2, i4),
                                  t2.indicesToPos({i1, i2, i4, i3}));
                    }
                }
            }
        }
    }


} // namespace Hammer::MultiDimensional
