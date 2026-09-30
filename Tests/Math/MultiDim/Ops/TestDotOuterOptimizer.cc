///
/// @file  TestDotOuterOptimizer.cc
/// @brief Tests for DotOuterOptimizer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"

#define DOTOUTEROPTIMIZER_FRIENDS             \
    FRIEND_TEST(DotOptimizerTest, OneAOneB);  \
    FRIEND_TEST(DotOptimizerTest, OneAOneC);  \
    FRIEND_TEST(DotOptimizerTest, OneATwoA);  \
    FRIEND_TEST(DotOptimizerTest, OneATwoB);  \
    FRIEND_TEST(DotOptimizerTest, OneATwoC);  \
    FRIEND_TEST(DotOptimizerTest, OneATwoD);  \
    FRIEND_TEST(DotOptimizerTest, OneAThree); \
    FRIEND_TEST(DotOptimizerTest, OneBOneC);  \
    FRIEND_TEST(DotOptimizerTest, OneBTwoA);  \
    FRIEND_TEST(DotOptimizerTest, OneBTwoB);  \
    FRIEND_TEST(DotOptimizerTest, OneBTwoC);  \
    FRIEND_TEST(DotOptimizerTest, OneBTwoD);  \
    FRIEND_TEST(DotOptimizerTest, OneBThree); \
    FRIEND_TEST(DotOptimizerTest, OneCTwoA);  \
    FRIEND_TEST(DotOptimizerTest, OneCTwoB);  \
    FRIEND_TEST(DotOptimizerTest, OneCTwoC);  \
    FRIEND_TEST(DotOptimizerTest, OneCTwoD);  \
    FRIEND_TEST(DotOptimizerTest, OneCThree); \
    FRIEND_TEST(DotOptimizerTest, TwoATwoB);  \
    FRIEND_TEST(DotOptimizerTest, TwoATwoC);  \
    FRIEND_TEST(DotOptimizerTest, TwoATwoD);  \
    FRIEND_TEST(DotOptimizerTest, TwoAThree); \
    FRIEND_TEST(DotOptimizerTest, TwoBTwoC);  \
    FRIEND_TEST(DotOptimizerTest, TwoBTwoD);  \
    FRIEND_TEST(DotOptimizerTest, TwoBThree); \
    FRIEND_TEST(DotOptimizerTest, TwoCTwoD);  \
    FRIEND_TEST(DotOptimizerTest, TwoCThree); \
    FRIEND_TEST(DotOptimizerTest, TwoDThree); \
    FRIEND_TEST(DotOptimizerTest, Other)

#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/Operations.hh"
#include "Hammer/Math/MultiDim/Ops/DotOuterOptimizer.hh"

#include "Hammer/Exceptions.hh"


using namespace std;

namespace Hammer::MultiDimensional::Ops {
    class DotOptimizerTest : public testing::Test {
    protected:

        // Base tensors available:
        //
        // A [WILSON_BCENU]
        // B [WILSON_BCENU_HC]
        // C [FF_BD]
        // D [FF_BD_HC]
        // E [FF_BD, WILSON_BCENU]
        // F [FF_BD_HC, WILSON_BCENU_HC]
        // G [WILSON_BCENU, FF_BD]
        // H [WILSON_BCENU_HC, FF_BD_HC]
        // I [FF_BD, FF_BD_VAR]
        // J [FF_BD_HC, FF_BD_VAR_HC]
        // K [FF_BD_VAR, FF_BD_VAR_HC]
        // L [WILSON_BCENU, WILSON_BCENU_HC]
        // M [WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC]
        //
        //
        // Combinations avialable
        //
        // 0) [[A, F], [B, F], [B, T], [A, T]]
        // 1) [[A, F], [A, T], [A, F], [A, T]]
        // 2) [[A, F], [L, T], [A, T]]
        //
        // 3) [[C, F], [L, F], [D, F]]
        // 4) [[C, F], [L, F], [C, T]]
        // 5) [[C, F], [A, F], [B, F], [D, F]]
        // 6) [[C, F], [A, F], [A, T], [C, T]]
        //
        // 7) [[M, F], [M, F]]
        //
        DotOptimizerTest() {
            auto v1a = makeVector({2}, {WILSON_BCENU}, {1., 2.i});
            auto v1ap = makeVector({2}, {WILSON_BCENU_HC}, {1., -2.i});
            auto v1b = makeVector({3}, {FF_BD}, {3.i, 4., 0.});
            auto v1bp = makeVector({3}, {FF_BD_HC}, {-3.i, 4., 0.});
            auto v2a = makeVector({3, 2}, {FF_BD, WILSON_BCENU}, {1., 0., 0., 0., 1.i, 1.});
            auto v2ap = makeVector({3, 2}, {FF_BD_HC, WILSON_BCENU_HC}, {1., 0., 0., 0., -1.i, 1.});
            auto v2b = makeVector({2, 3}, {WILSON_BCENU, FF_BD}, {1., 0., 0., 0., 1.i, 1.});
            auto v2bp = makeVector({2, 3}, {WILSON_BCENU_HC, FF_BD_HC}, {1., 0., 0., 0., -1.i, 1.});
            auto v2c = makeVector({3, 2}, {FF_BD, FF_BD_VAR}, {1., 1., 1., 2.i, 2.i, 2.i});
            auto v2cp = makeVector({3, 2}, {FF_BD_HC, FF_BD_VAR_HC}, {1., 1., 1., -2.i, -2.i, -2.i});
            auto v2d = makeVector({2, 2}, {FF_BD_VAR, FF_BD_VAR_HC}, {1., 0, 0, 1.});
            auto v2e = makeVector({2, 2}, {WILSON_BCENU, WILSON_BCENU_HC}, {1., 0, 0, 1.});
            auto v3 = makeVector({2, 3, 2, 3}, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC},
                                 {1.,  0, 0, 1.,  2., -3., 5.,  0, 0, 5.,  2.i, -3.i, 7.,  0, 0, 5.,  2.i, -3.,
                                  11., 0, 0, 11., 2., -3., 13., 0, 0, 13., 2.i, -3.i, 17., 0, 0, 17., 2.i, -3.});
            _data.emplace_back(SharedTensorData{v1a.release()});
            _data.emplace_back(SharedTensorData{v1ap.release()});
            _data.emplace_back(SharedTensorData{v1b.release()});
            _data.emplace_back(SharedTensorData{v1bp.release()});
            _data.emplace_back(SharedTensorData{v2a.release()});
            _data.emplace_back(SharedTensorData{v2ap.release()});
            _data.emplace_back(SharedTensorData{v2b.release()});
            _data.emplace_back(SharedTensorData{v2bp.release()});
            _data.emplace_back(SharedTensorData{v2c.release()});
            _data.emplace_back(SharedTensorData{v2cp.release()});
            _data.emplace_back(SharedTensorData{v2d.release()});
            _data.emplace_back(SharedTensorData{v2e.release()});
            _data.emplace_back(SharedTensorData{v3.release()});

            vector<pair<SharedTensorData, bool>> vecL1a; // {WC, WC_HC, WC, WC_HC} two tensors
            vecL1a.emplace_back(_data[0], false);
            vecL1a.emplace_back(_data[1], false);
            vecL1a.emplace_back(_data[1], true);
            vecL1a.emplace_back(_data[0], true);
            BlockIndexing b1a{{{2}, {2}, {2}, {2}},
                              {{WILSON_BCENU}, {WILSON_BCENU_HC}, {WILSON_BCENU}, {WILSON_BCENU_HC}}};
            auto o1a = combineSharedTensors(std::move(vecL1a));
            _tensors.emplace_back(std::move(o1a));
            _indices.emplace_back(std::move(b1a));

            vector<pair<SharedTensorData, bool>> vecL1b; // {WC, WC_HC, WC, WC_HC} one tensor
            vecL1b.emplace_back(_data[0], false);
            vecL1b.emplace_back(_data[0], true);
            vecL1b.emplace_back(_data[0], false);
            vecL1b.emplace_back(_data[0], true);
            BlockIndexing b1b{{{2}, {2}, {2}, {2}},
                              {{WILSON_BCENU}, {WILSON_BCENU_HC}, {WILSON_BCENU}, {WILSON_BCENU_HC}}};
            auto o1b = combineSharedTensors(std::move(vecL1b));
            _tensors.emplace_back(std::move(o1b));
            _indices.emplace_back(std::move(b1b));

            vector<pair<SharedTensorData, bool>> vecL1c; // {WC, [WC, WC_HC], WC_HC} two tensors
            vecL1c.emplace_back(_data[0], false);
            vecL1c.emplace_back(_data[11], true);
            vecL1c.emplace_back(_data[0], true);
            BlockIndexing b1c{{{2}, {2, 2}, {2}}, {{WILSON_BCENU}, {WILSON_BCENU_HC, WILSON_BCENU}, {WILSON_BCENU_HC}}};
            auto o1c = combineSharedTensors(std::move(vecL1c));
            _tensors.emplace_back(std::move(o1c));
            _indices.emplace_back(std::move(b1c));

            vector<pair<SharedTensorData, bool>> vecL2a; // {FF, [WC, WC_HC], FF_HC} three tensors
            vecL2a.emplace_back(_data[2], false);
            vecL2a.emplace_back(_data[11], false);
            vecL2a.emplace_back(_data[3], false);
            auto o2a = combineSharedTensors(std::move(vecL2a));
            BlockIndexing b2a{{{3}, {2, 2}, {3}}, {{FF_BD}, {WILSON_BCENU, WILSON_BCENU_HC}, {FF_BD_HC}}};
            _tensors.emplace_back(std::move(o2a));
            _indices.emplace_back(std::move(b2a));

            vector<pair<SharedTensorData, bool>> vecL2b; // {FF, [WC, WC_HC], FF_HC} two tensors
            vecL2b.emplace_back(_data[2], false);
            vecL2b.emplace_back(_data[11], false);
            vecL2b.emplace_back(_data[2], true);
            BlockIndexing b2b{{{3}, {2, 2}, {3}}, {{FF_BD}, {WILSON_BCENU, WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2b = combineSharedTensors(std::move(vecL2b));
            _tensors.emplace_back(std::move(o2b));
            _indices.emplace_back(std::move(b2b));

            vector<pair<SharedTensorData, bool>> vecL2c; // {FF, WC, WC_HC, FF_HC} 4 tensors
            vecL2c.emplace_back(_data[2], false);
            vecL2c.emplace_back(_data[0], false);
            vecL2c.emplace_back(_data[1], false);
            vecL2c.emplace_back(_data[3], false);
            BlockIndexing b2c{{{3}, {2}, {2}, {3}}, {{FF_BD}, {WILSON_BCENU}, {WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2c = combineSharedTensors(std::move(vecL2c));
            _tensors.emplace_back(std::move(o2c));
            _indices.emplace_back(std::move(b2c));

            vector<pair<SharedTensorData, bool>> vecL2d; // {FF, WC, WC_HC, FF_HC} 2 tensors
            vecL2d.emplace_back(_data[2], false);
            vecL2d.emplace_back(_data[0], false);
            vecL2d.emplace_back(_data[0], true);
            vecL2d.emplace_back(_data[2], true);
            BlockIndexing b2d{{{3}, {2}, {2}, {3}}, {{FF_BD}, {WILSON_BCENU}, {WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2d = combineSharedTensors(std::move(vecL2d));
            _tensors.emplace_back(std::move(o2d));
            _indices.emplace_back(std::move(b2d));

            vector<pair<SharedTensorData, bool>> vecL3; // {WC, FF, WC_HC, FF_HC, WC, FF, WC_HC, FF_HC} 2 tensors
            vecL3.emplace_back(_data[12], false);
            vecL3.emplace_back(_data[12], false);
            BlockIndexing b3{
                {{2, 3, 2, 3}, {2, 3, 2, 3}},
                {{WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC}, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC}}};
            auto o3 = combineSharedTensors(std::move(vecL3));
            _tensors.emplace_back(std::move(o3));
            _indices.emplace_back(std::move(b3));
        }


        std::vector<SharedTensorData> _data;
        std::vector<TensorData> _tensors;
        std::vector<BlockIndexing> _indices;
    };

    // static const IContainer* getIC(SharedTensorData& p) {
    //     return const_cast<const IContainer*>(p.get());
    // }

    // static const IContainer* getIC(TensorData& p) {
    //     return const_cast<const IContainer*>(p.get());
    // }

    // static OuterContainer* getOC(SharedTensorData& o) {
    //     return static_cast<OuterContainer*>(o.get());
    // }

    static OuterContainer* getOC(TensorData& o) {
        return static_cast<OuterContainer*>(o.get());
    }

    TEST_F(DotOptimizerTest, OneAOneB) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[1]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 8);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 4);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        {
            auto it = findEdge(2, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        {
            auto it = findEdge(3, 64 + 3);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }


    TEST_F(DotOptimizerTest, OneAOneC) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[2]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 3);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        {
            auto it = findEdge(3, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._stars[0]._isCenterLeft, false);
        EXPECT_EQ(optim._stars[0]._edges.size(), 2);
        // After sorted edge order: edges (0,64),(1,65),(2,65),(3,66)
        // Star on RIGHT 65: EdgeId=1 (from=1) and EdgeId=2 (from=2)
        // EdgeId=1 processed first -> EquivId=0; EdgeId=2 -> EquivId=1
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[0].first]._fromVertex, 1);
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[0].first]._toVertex, 64 + 1);
        EXPECT_EQ(optim._stars[0]._edges[0].second, 0);
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[1].first]._fromVertex, 2);
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[1].first]._toVertex, 64 + 1);
        EXPECT_EQ(optim._stars[0]._edges[1].second, 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneATwoA) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[3]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 5);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._stars[0]._isCenterLeft, false);
        EXPECT_EQ(optim._stars[0]._edges.size(), 2);
        // After sorted edge order: edges (0,65),(1,65)
        // Star on RIGHT 65: EdgeId=0 (from=0, WC pos 0) -> EquivId=0; EdgeId=1 (from=1, WC_HC pos 1) ->
        // EquivId=1
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[0].first]._fromVertex, 0);
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[0].first]._toVertex, 64 + 1);
        EXPECT_EQ(optim._stars[0]._edges[0].second, 0);
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[1].first]._fromVertex, 1);
        EXPECT_EQ(optim._edges[optim._stars[0]._edges[1].first]._toVertex, 64 + 1);
        EXPECT_EQ(optim._stars[0]._edges[1].second, 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, OneATwoB) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 5);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, OneATwoC) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 8);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 6);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneATwoD) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 8);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 6);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneAThree) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 6);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 2);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        {
            auto it = findEdge(2, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(3, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 2);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneBOneC) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[2]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 3);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        {
            auto it = findEdge(3, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneBTwoA) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[3]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 5);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, OneBTwoB) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 5);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, OneBTwoC) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 8);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 6);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneBTwoD) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 8);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 6);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, OneBThree) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 6);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 2);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        {
            auto it = findEdge(2, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(3, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 2);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, OneCTwoA) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[3]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 6);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 4);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneCTwoB) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 6);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 4);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneCTwoC) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 5);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneCTwoD) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 2);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 5);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        EXPECT_EQ(optim._vertices[2]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 0]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._vertices[64 + 3]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, OneCThree) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 5);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 1);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 2);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 4);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 4);
    }

    TEST_F(DotOptimizerTest, TwoATwoB) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 3);
        EXPECT_EQ(optim._vertices.size(), 6);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 3);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 2);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
            EXPECT_EQ(it->_localPositions.back().first, 1);
            EXPECT_EQ(it->_localPositions.back().second, 1);
        }
        {
            auto it = findEdge(2, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 3);
    }

    TEST_F(DotOptimizerTest, TwoATwoC) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 3);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 3);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 3);
    }

    TEST_F(DotOptimizerTest, TwoATwoD) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 3);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 3);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 3);
    }

    TEST_F(DotOptimizerTest, TwoAThree) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 3);
        EXPECT_EQ(optim._vertices.size(), 5);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 2);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        {
            auto it = findEdge(1, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 2);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
            EXPECT_EQ(it->_localPositions.back().first, 1);
            EXPECT_EQ(it->_localPositions.back().second, 2);
        }
        {
            auto it = findEdge(2, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 3);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._vertices[64 + 1]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, TwoBTwoC) {
        auto* one = getOC(_tensors[4]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 3);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 3);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 3);
    }

    TEST_F(DotOptimizerTest, TwoBTwoD) {
        auto* one = getOC(_tensors[4]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 7);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 3);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().first, 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 3);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, TwoBThree) {
        auto* one = getOC(_tensors[4]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 3);
        EXPECT_EQ(optim._vertices.size(), 5);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 2);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        {
            auto it = findEdge(1, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 2);
            EXPECT_EQ(it->_localPositions.front().first, 0);
            EXPECT_EQ(it->_localPositions.front().second, 0);
            EXPECT_EQ(it->_localPositions.back().first, 1);
            EXPECT_EQ(it->_localPositions.back().second, 2);
        }
        {
            auto it = findEdge(2, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 3);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._vertices[64 + 1]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 2);
    }

    TEST_F(DotOptimizerTest, TwoCTwoD) {
        auto* one = getOC(_tensors[5]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 8);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 4);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64 + 0);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(1, 64 + 1);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64 + 2);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(3, 64 + 3);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 2);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 4);
    }

    TEST_F(DotOptimizerTest, TwoCThree) {
        auto* one = getOC(_tensors[5]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 6);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 2);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        {
            auto it = findEdge(1, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        {
            auto it = findEdge(3, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 3);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        EXPECT_EQ(optim._vertices[64 + 1]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, TwoDThree) {
        auto* one = getOC(_tensors[6]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto optim = Ops::DotOuterOptimizer(*one, *other, idxs);
        EXPECT_EQ(optim._edges.size(), 4);
        EXPECT_EQ(optim._vertices.size(), 6);
        EXPECT_EQ(optim._stars.size(), 0);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 0);
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 0);
        EXPECT_EQ(optim._maxTensorCount, 2);
        EXPECT_EQ(optim._nextVertexId, 0);
        auto findEdge = [&](DotOuterOptimizer::VertexId from, DotOuterOptimizer::VertexId to) {
            return std::find_if(optim._edges.begin(), optim._edges.end(),
                                [from, to](const auto& e) { return e._fromVertex == from && e._toVertex == to; });
        };
        {
            auto it = findEdge(0, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 1);
        }
        {
            auto it = findEdge(1, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 0);
        }
        {
            auto it = findEdge(2, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 2);
        }
        {
            auto it = findEdge(3, 64);
            ASSERT_NE(it, optim._edges.end());
            EXPECT_EQ(it->_localPositions.size(), 1);
            EXPECT_EQ(it->_localPositions.front().second, 3);
        }
        optim.fillEquivalentSubGrGroups();
        EXPECT_EQ(optim._stars.size(), 1);
        EXPECT_EQ(optim._equivalentSubGrGroups.size(), 1);
        EXPECT_EQ(optim._vertices[64 + 1]._type, DotOuterOptimizer::Vertex::Type::UNUSED);
        optim.assignDataToVertices(*(one->begin()), *(other->begin()));
        optim.createIdenticalSubGrGroups();
        EXPECT_EQ(optim._identicalSubGrGroups.size(), 1);
    }

    TEST_F(DotOptimizerTest, Other) {
        auto v1 = makeVector({2}, {WILSON_BCENU}, {1.i, 2.});
        auto v2 = makeVector({2}, {WILSON_BCENU_HC}, {-2., 1.});
        auto v3 = makeVector({3, 2}, {FF_BD, WILSON_BCENU}, {-2., 1., 1.i, 5., 3., 1.});
        auto v4 = makeVector({3, 2}, {FF_BD_HC, WILSON_BCENU_HC}, {3., -2.i, 4., 1., -2., 1.});

        vector<TensorData> vec1;
        vec1.push_back(std::move(v1));
        vec1.push_back(std::move(v2));
        vector<TensorData> vec2;
        vec2.push_back(std::move(v3));
        vec2.push_back(std::move(v4));

        auto o1 = combineTensors(vec1);
        auto o2 = combineTensors(vec2);

        auto idxs2 = o1->getSameLabelPairs(*o2, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});

        auto optim2 = Ops::DotOuterOptimizer(*getOC(o1), *getOC(o2), idxs2);
        optim2.fillEquivalentSubGrGroups();
        optim2.assignDataToVertices(*(getOC(o1)->begin()), *(getOC(o2)->begin()));
        optim2.createIdenticalSubGrGroups();
        EXPECT_EQ(optim2._stars.size(), 2);
        EXPECT_EQ(optim2._equivalentSubGrGroups.size(), 1);
        EXPECT_EQ(optim2._identicalSubGrGroups.size(), 2);
    }

} // namespace Hammer::MultiDimensional::Ops
