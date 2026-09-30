///
/// @file  TestDotOuter.cc
/// @brief Tests for DotOuter
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/Operations.hh"
#include "Hammer/Math/MultiDim/Ops/Dot.hh"
#include "Hammer/Math/MultiDim/BruteForceIterator.hh"
#include "Hammer/Exceptions.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    class DotOuterTest : public testing::Test {
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
        DotOuterTest() {
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

            for (auto& elem : _data) {
                _dataS.emplace_back(SharedTensorData{toSparse(elem->clone()).release()});
            }

            vector<pair<SharedTensorData, bool>> vecL1a; // {WC, WC_HC, WC_HC, WC} two tensors
            vecL1a.emplace_back(_data[0], false);
            vecL1a.emplace_back(_data[1], false);
            vecL1a.emplace_back(_data[1], true);
            vecL1a.emplace_back(_data[0], true);
            vector<pair<SharedTensorData, bool>> sparL1a; // {WC, WC_HC, WC_HC, WC} two tensors
            sparL1a.emplace_back(_dataS[0], false);
            sparL1a.emplace_back(_dataS[1], false);
            sparL1a.emplace_back(_dataS[1], true);
            sparL1a.emplace_back(_dataS[0], true);
            BlockIndexing b1a{{{2}, {2}, {2}, {2}},
                              {{WILSON_BCENU}, {WILSON_BCENU_HC}, {WILSON_BCENU}, {WILSON_BCENU_HC}}};
            auto o1a = combineSharedTensors(std::move(vecL1a));
            _tensors.emplace_back(std::move(o1a));
            auto os1a = combineSharedTensors(std::move(sparL1a));
            _tensorsS.emplace_back(std::move(os1a));
            _indices.emplace_back(std::move(b1a));

            vector<pair<SharedTensorData, bool>> vecL1b; // {WC, WC_HC, WC_HC, WC} one tensor
            vecL1b.emplace_back(_data[0], false);
            vecL1b.emplace_back(_data[0], true);
            vecL1b.emplace_back(_data[0], false);
            vecL1b.emplace_back(_data[0], true);
            vector<pair<SharedTensorData, bool>> sparL1b; // {WC, WC_HC, WC_HC, WC} one tensor
            sparL1b.emplace_back(_dataS[0], false);
            sparL1b.emplace_back(_dataS[0], true);
            sparL1b.emplace_back(_dataS[0], false);
            sparL1b.emplace_back(_dataS[0], true);
            BlockIndexing b1b{{{2}, {2}, {2}, {2}},
                              {{WILSON_BCENU}, {WILSON_BCENU_HC}, {WILSON_BCENU}, {WILSON_BCENU_HC}}};
            auto o1b = combineSharedTensors(std::move(vecL1b));
            _tensors.emplace_back(std::move(o1b));
            auto os1b = combineSharedTensors(std::move(sparL1b));
            _tensorsS.emplace_back(std::move(os1b));
            _indices.emplace_back(std::move(b1b));

            vector<pair<SharedTensorData, bool>> vecL1c; // {WC, [WC, WC_HC], WC_HC} two tensors
            vecL1c.emplace_back(_data[0], false);
            vecL1c.emplace_back(_data[11], true);
            vecL1c.emplace_back(_data[0], true);
            vector<pair<SharedTensorData, bool>> sparL1c; // {WC, [WC, WC_HC], WC_HC} two tensors
            sparL1c.emplace_back(_dataS[0], false);
            sparL1c.emplace_back(_dataS[11], true);
            sparL1c.emplace_back(_dataS[0], true);
            BlockIndexing b1c{{{2}, {2, 2}, {2}}, {{WILSON_BCENU}, {WILSON_BCENU_HC, WILSON_BCENU}, {WILSON_BCENU_HC}}};
            auto o1c = combineSharedTensors(std::move(vecL1c));
            _tensors.emplace_back(std::move(o1c));
            auto os1c = combineSharedTensors(std::move(sparL1c));
            _tensorsS.emplace_back(std::move(os1c));
            _indices.emplace_back(std::move(b1c));

            vector<pair<SharedTensorData, bool>> vecL2a; // {FF, [WC, WC_HC], FF_HC} three tensors
            vecL2a.emplace_back(_data[2], false);
            vecL2a.emplace_back(_data[11], false);
            vecL2a.emplace_back(_data[3], false);
            vector<pair<SharedTensorData, bool>> sparL2a; // {FF, [WC, WC_HC], FF_HC} three tensors
            sparL2a.emplace_back(_dataS[2], false);
            sparL2a.emplace_back(_dataS[11], false);
            sparL2a.emplace_back(_dataS[3], false);
            BlockIndexing b2a{{{3}, {2, 2}, {3}}, {{FF_BD}, {WILSON_BCENU, WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2a = combineSharedTensors(std::move(vecL2a));
            _tensors.emplace_back(std::move(o2a));
            auto os2a = combineSharedTensors(std::move(sparL2a));
            _tensorsS.emplace_back(std::move(os2a));
            _indices.emplace_back(std::move(b2a));

            vector<pair<SharedTensorData, bool>> vecL2b; // {FF, [WC, WC_HC], FF_HC} two tensors
            vecL2b.emplace_back(_data[2], false);
            vecL2b.emplace_back(_data[11], false);
            vecL2b.emplace_back(_data[2], true);
            vector<pair<SharedTensorData, bool>> sparL2b; // {FF, [WC, WC_HC], FF_HC} two tensors
            sparL2b.emplace_back(_dataS[2], false);
            sparL2b.emplace_back(_dataS[11], false);
            sparL2b.emplace_back(_dataS[2], true);
            BlockIndexing b2b{{{3}, {2, 2}, {3}}, {{FF_BD}, {WILSON_BCENU, WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2b = combineSharedTensors(std::move(vecL2b));
            _tensors.emplace_back(std::move(o2b));
            auto os2b = combineSharedTensors(std::move(sparL2b));
            _tensorsS.emplace_back(std::move(os2b));
            _indices.emplace_back(std::move(b2b));

            vector<pair<SharedTensorData, bool>> vecL2c; // {FF, WC, WC_HC, FF_HC} 4 tensors
            vecL2c.emplace_back(_data[2], false);
            vecL2c.emplace_back(_data[0], false);
            vecL2c.emplace_back(_data[1], false);
            vecL2c.emplace_back(_data[3], false);
            vector<pair<SharedTensorData, bool>> sparL2c; // {FF, WC, WC_HC, FF_HC} 4 tensors
            sparL2c.emplace_back(_dataS[2], false);
            sparL2c.emplace_back(_dataS[0], false);
            sparL2c.emplace_back(_dataS[1], false);
            sparL2c.emplace_back(_dataS[3], false);
            BlockIndexing b2c{{{3}, {2}, {2}, {3}}, {{FF_BD}, {WILSON_BCENU}, {WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2c = combineSharedTensors(std::move(vecL2c));
            _tensors.emplace_back(std::move(o2c));
            auto os2c = combineSharedTensors(std::move(sparL2c));
            _tensorsS.emplace_back(std::move(os2c));
            _indices.emplace_back(std::move(b2c));

            vector<pair<SharedTensorData, bool>> vecL2d; // {FF, WC, WC_HC, FF_HC} 2 tensors
            vecL2d.emplace_back(_data[2], false);
            vecL2d.emplace_back(_data[0], false);
            vecL2d.emplace_back(_data[0], true);
            vecL2d.emplace_back(_data[2], true);
            vector<pair<SharedTensorData, bool>> sparL2d; // {FF, WC, WC_HC, FF_HC} 2 tensors
            sparL2d.emplace_back(_dataS[2], false);
            sparL2d.emplace_back(_dataS[0], false);
            sparL2d.emplace_back(_dataS[0], true);
            sparL2d.emplace_back(_dataS[2], true);
            BlockIndexing b2d{{{3}, {2}, {2}, {3}}, {{FF_BD}, {WILSON_BCENU}, {WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2d = combineSharedTensors(std::move(vecL2d));
            _tensors.emplace_back(std::move(o2d));
            auto os2d = combineSharedTensors(std::move(sparL2d));
            _tensorsS.emplace_back(std::move(os2d));
            _indices.emplace_back(std::move(b2d));

            vector<pair<SharedTensorData, bool>> vecL3; // {WC, FF, WC_HC, FF_HC, WC, FF, WC_HC, FF_HC} 2 tensors
            vecL3.emplace_back(_data[12], false);
            vecL3.emplace_back(_data[12], false);
            vector<pair<SharedTensorData, bool>> sparL3; // {WC, FF, WC_HC, FF_HC, WC, FF, WC_HC, FF_HC} 2 tensors
            sparL3.emplace_back(_dataS[12], false);
            sparL3.emplace_back(_dataS[12], false);
            BlockIndexing b3{
                {{2, 3, 2, 3}, {2, 3, 2, 3}},
                {{WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC}, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC}}};
            auto o3 = combineSharedTensors(std::move(vecL3));
            _tensors.emplace_back(std::move(o3));
            auto os3 = combineSharedTensors(std::move(sparL3));
            _tensorsS.emplace_back(std::move(os3));
            _indices.emplace_back(std::move(b3));
        }


        std::vector<SharedTensorData> _data;
        std::vector<SharedTensorData> _dataS;
        std::vector<TensorData> _tensors;
        std::vector<TensorData> _tensorsS;
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

    TEST_F(DotOuterTest, ContractSinglesVV) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[1]);
        auto res = TensorData{Ops::Dot::contractSingles(*(one->begin()->begin()), *(other->begin()->begin()))};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractSingles(one->begin()->at(2), *(other->begin()->begin()))};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }

    TEST_F(DotOuterTest, ContractSinglesSS) {
        auto* one = getOC(_tensorsS[0]);
        auto* other = getOC(_tensorsS[1]);
        auto res = TensorData{Ops::Dot::contractSingles(*(one->begin()->begin()), *(other->begin()->begin()))};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractSingles(one->begin()->at(2), *(other->begin()->begin()))};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }

    TEST_F(DotOuterTest, ContractSinglesSV) {
        auto* one = getOC(_tensorsS[0]);
        auto* other = getOC(_tensors[1]);
        auto res = TensorData{Ops::Dot::contractSingles(*(one->begin()->begin()), *(other->begin()->begin()))};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractSingles(one->begin()->at(2), *(other->begin()->begin()))};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }

    TEST_F(DotOuterTest, ContractSinglesVS) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensorsS[1]);
        auto res = TensorData{Ops::Dot::contractSingles(*(one->begin()->begin()), *(other->begin()->begin()))};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractSingles(one->begin()->at(2), *(other->begin()->begin()))};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), -3.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }

    TEST_F(DotOuterTest, ContractStarFullVV) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother2, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }


    TEST_F(DotOuterTest, ContractStarFullSV) {
        auto* one = getOC(_tensorsS[3]);
        auto* other = getOC(_tensors[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother2, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }


    TEST_F(DotOuterTest, ContractStarFullVS) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensorsS[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother2, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }


    TEST_F(DotOuterTest, ContractStarFullSS) {
        auto* one = getOC(_tensorsS[3]);
        auto* other = getOC(_tensorsS[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res->rank(), 0);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res->element({}).imag(), 0.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(1), vecother2, {{0, 0}, {1, 1}})};
        EXPECT_EQ(res2->rank(), 0);
        EXPECT_DOUBLE_EQ(res2->element({}).real(), 5.);
        EXPECT_DOUBLE_EQ(res2->element({}).imag(), 0.);
    }


    TEST_F(DotOuterTest, ContractStarPartialVV) {
        auto* one = getOC(_tensors[7]);
        auto* other = getOC(_tensors[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res->rank(), 2);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 45.);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 20.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).real(), 8.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).imag(), -4.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).imag(), 6.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 57.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 16.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).real(), -6.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).imag(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 75.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 24.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).imag(), 6.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother2, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res2->rank(), 2);
        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                EXPECT_DOUBLE_EQ(res2->element({i, j}).real(), res->element({i, j}).real());
                EXPECT_DOUBLE_EQ(res2->element({i, j}).imag(), res->element({i, j}).imag());
            }
        }
    }


    TEST_F(DotOuterTest, ContractStarPartialVS) {
        auto* one = getOC(_tensors[7]);
        auto* other = getOC(_tensorsS[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res->rank(), 2);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 45.);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 20.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).real(), 8.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).imag(), -4.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).imag(), 6.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 57.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 16.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).real(), -6.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).imag(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 75.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 24.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).imag(), 6.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother2, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res2->rank(), 2);
        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                EXPECT_DOUBLE_EQ(res2->element({i, j}).real(), res->element({i, j}).real());
                EXPECT_DOUBLE_EQ(res2->element({i, j}).imag(), res->element({i, j}).imag());
            }
        }
    }


    TEST_F(DotOuterTest, ContractStarPartialSV) {
        auto* one = getOC(_tensorsS[7]);
        auto* other = getOC(_tensors[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res->rank(), 2);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 45.);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 20.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).real(), 8.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).imag(), -4.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).imag(), 6.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 57.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 16.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).real(), -6.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).imag(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 75.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 24.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).imag(), 6.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother2, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res2->rank(), 2);
        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                EXPECT_DOUBLE_EQ(res2->element({i, j}).real(), res->element({i, j}).real());
                EXPECT_DOUBLE_EQ(res2->element({i, j}).imag(), res->element({i, j}).imag());
            }
        }
    }


    TEST_F(DotOuterTest, ContractStarPartialSS) {
        auto* one = getOC(_tensorsS[7]);
        auto* other = getOC(_tensorsS[0]);
        auto vecother = {other->begin()->at(0), other->begin()->at(3)};
        auto vecother2 = {other->begin()->at(0), other->begin()->at(1)};
        auto res = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res->rank(), 2);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 45.);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 20.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).real(), 8.);
        EXPECT_DOUBLE_EQ(res->element({0, 1}).imag(), -4.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({0, 2}).imag(), 6.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 57.);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 16.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({1, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).real(), -6.);
        EXPECT_DOUBLE_EQ(res->element({1, 2}).imag(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 75.);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 24.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(res->element({2, 1}).imag(), 8.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).real(), -12.);
        EXPECT_DOUBLE_EQ(res->element({2, 2}).imag(), 6.);
        auto res2 = TensorData{Ops::Dot::contractStar(one->begin()->at(0), vecother2, {{0, 0}, {2, 1}})};
        EXPECT_EQ(res2->rank(), 2);
        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                EXPECT_DOUBLE_EQ(res2->element({i, j}).real(), res->element({i, j}).real());
                EXPECT_DOUBLE_EQ(res2->element({i, j}).imag(), res->element({i, j}).imag());
            }
        }
    }


    // Helper: run Dot(a, b) through the direct VectorContainer×OuterContainer overload
    static TensorData runDotVO(VectorContainer& a, OuterContainer& b, const IndexPairList& idxs,
                               pair<bool, bool> hc = {false, false}) {
        Ops::Dot dotter{idxs, hc};
        return TensorData{dotter(a, b)};
    }

    // Helper: run the same Dot through the BruteForce (IContainer, IContainer) fallback
    static TensorData runBrute(IContainer& a, IContainer& b, const IndexPairList& idxs,
                               pair<bool, bool> hc = {false, false}) {
        Ops::Dot dotter{idxs, hc};
        return TensorData{dotter.operator()(a, static_cast<const IContainer&>(b))};
    }

    static void expectTensorEq(const IContainer& got, const IContainer& ref) {
        ASSERT_EQ(got.rank(), ref.rank());
        BruteForceIterator bf{ref.dims()};
        for (auto idx : bf) {
            EXPECT_DOUBLE_EQ(got.element(idx).real(), ref.element(idx).real());
            EXPECT_DOUBLE_EQ(got.element(idx).imag(), ref.element(idx).imag());
        }
    }

    TEST_F(DotOuterTest, VectorOuterStarPartial) {
        // [E, F] . [[C, F], [A, F], [B, F], [D, F]]
        auto* a = static_cast<VectorContainer*>(_data[4].get());
        auto* b = getOC(_tensors[5]);
        auto idxs = a->getSameLabelPairs(*b, {WILSON_BCENU});
        auto fast = runDotVO(*a, *b, idxs);
        auto ref = runBrute(*a, *b, idxs);
        ASSERT_EQ(fast->rank(), ref->rank());
        expectTensorEq(*fast, *ref);
    }

    TEST_F(DotOuterTest, VectorOuterStarScalar) {
        // [L, F] . [[A, F], [B, F]]
        auto* a = static_cast<VectorContainer*>(_data[11].get());
        vector<pair<SharedTensorData, bool>> facs;
        facs.emplace_back(_data[0], false);
        facs.emplace_back(_data[1], false);
        auto bData = combineSharedTensors(std::move(facs));
        auto* b = getOC(bData);
        auto idxs = a->getSameLabelPairs(*b, {WILSON_BCENU, WILSON_BCENU_HC});
        auto fast = runDotVO(*a, *b, idxs);
        auto ref = runBrute(*a, *b, idxs);
        ASSERT_EQ(fast->rank(), 0u);
        EXPECT_DOUBLE_EQ(fast->element({}).real(), ref->element({}).real());
        EXPECT_DOUBLE_EQ(fast->element({}).imag(), ref->element({}).imag());
    }

    TEST_F(DotOuterTest, VectorOuterBoomerang) {
        // [L,F] . [[G,F],[H,F]]
        auto* a = static_cast<VectorContainer*>(_data[11].get());
        vector<pair<SharedTensorData, bool>> facs;
        facs.emplace_back(_data[6], false);
        facs.emplace_back(_data[7], false);
        auto bData = combineSharedTensors(std::move(facs));
        auto* b = getOC(bData);
        auto idxs = a->getSameLabelPairs(*b, {WILSON_BCENU, WILSON_BCENU_HC});
        auto fast = runDotVO(*a, *b, idxs);
        auto ref = runBrute(*a, *b, idxs);
        ASSERT_EQ(fast->rank(), ref->rank());
        expectTensorEq(*fast, *ref);
    }

    TEST_F(DotOuterTest, VectorOuterGeneral) {
        // [E,F] . [[C, F], [L, F], [D, F]]
        auto* a = static_cast<VectorContainer*>(_data[4].get());
        auto* b = getOC(_tensors[3]);
        auto idxs = a->getSameLabelPairs(*b, {FF_BD, WILSON_BCENU});
        auto fast = runDotVO(*a, *b, idxs);
        auto ref = runBrute(*a, *b, idxs);
        ASSERT_EQ(fast->rank(), ref->rank());
        expectTensorEq(*fast, *ref);
    }

    TEST_F(DotOuterTest, VectorOuterScalar) {
        // [M,F] . [[A, F], [C, F], [B, F], [D, F]]
        auto* a = static_cast<VectorContainer*>(_data[12].get());
        vector<pair<SharedTensorData, bool>> vec;
        vec.emplace_back(_data[0], false);
        vec.emplace_back(_data[2], false);
        vec.emplace_back(_data[1], false);
        vec.emplace_back(_data[3], false);
        auto data = combineSharedTensors(std::move(vec));
        auto* b = getOC(data);
        auto idxs = a->getSameLabelPairs(*b, {WILSON_BCENU, FF_BD, WILSON_BCENU_HC, FF_BD_HC});
        auto fast = runDotVO(*a, *b, idxs);
        auto ref = runBrute(*a, *b, idxs);
        ASSERT_EQ(fast->rank(), 0u);
        EXPECT_DOUBLE_EQ(fast->element({}).real(), ref->element({}).real());
        EXPECT_DOUBLE_EQ(fast->element({}).imag(), ref->element({}).imag());
    }

} // namespace Hammer::MultiDimensional
