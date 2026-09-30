///
/// @file  TestDotOuter2.cc
/// @brief Tests for DotOuter2
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
        // 1a) [[A, F], [B, F], [B, T], [A, T]]
        // 1b) [[A, F], [A, T], [A, F], [A, T]]
        // 1c) [[A, F], [L, T], [A, T]]
        //
        // 2a) [[C, F], [L, F], [D, F]]
        // 2b) [[C, F], [L, F], [C, T]]
        // 2c) [[C, F], [A, F], [B, F], [D, F]]
        // 2d) [[C, F], [A, F], [A, T], [C, T]]
        //
        // 3) [[M, F], [M, F]]
        //
        //
        DotOuterTest() {
            auto v1a = makeVector({2}, {WILSON_BCENU}, {1., 2.i});
            auto v1ap = makeVector({2}, {WILSON_BCENU_HC}, {1., -3.i});
            auto v1b = makeVector({3}, {FF_BD}, {3.i, 4., 0.});
            auto v1bp = makeVector({3}, {FF_BD_HC}, {-3.i, 5., 0.});
            auto v2a = makeVector({3, 2}, {FF_BD, WILSON_BCENU}, {1., 0., 0., 0., 1.i, 1.});
            auto v2ap = makeVector({3, 2}, {FF_BD_HC, WILSON_BCENU_HC}, {1., 0., 0., 0., -1.i, 2.});
            auto v2b = makeVector({2, 3}, {WILSON_BCENU, FF_BD}, {1., 0., 2., 0., 1.i, 1.});
            auto v2bp = makeVector({2, 3}, {WILSON_BCENU_HC, FF_BD_HC}, {3., 0., 0., 0., -1.i, 1.});
            auto v2c = makeVector({3, 2}, {FF_BD, FF_BD_VAR}, {1., 1., 1., 2.i, 2.i, 2.i});
            auto v2cp = makeVector({3, 2}, {FF_BD_HC, FF_BD_VAR_HC}, {1., 5., 1., -2.i, -2.i, -4.i});
            auto v2d = makeVector({2, 2}, {FF_BD_VAR, FF_BD_VAR_HC}, {1., 0, 0, 1.});
            auto v2e = makeVector({2, 2}, {WILSON_BCENU, WILSON_BCENU_HC}, {3., 0, 0, 3.});
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

            vector<pair<SharedTensorData, bool>> vecL1a; // {WC, WC_HC, WC_HC, WC} two tensors
            vecL1a.emplace_back(_data[0], false);
            vecL1a.emplace_back(_data[1], false);
            vecL1a.emplace_back(_data[1], true);
            vecL1a.emplace_back(_data[0], true);
            BlockIndexing b1a{{{2}, {2}, {2}, {2}},
                              {{WILSON_BCENU}, {WILSON_BCENU_HC}, {WILSON_BCENU}, {WILSON_BCENU_HC}}};
            auto o1a = combineSharedTensors(std::move(vecL1a));
            _tensors.emplace_back(std::move(o1a));
            _indices.emplace_back(std::move(b1a));

            vector<pair<SharedTensorData, bool>> vecL1b; // {WC, WC_HC, WC_HC, WC} one tensor
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
            BlockIndexing b2a{{{3}, {2, 2}, {3}}, {{FF_BD}, {WILSON_BCENU, WILSON_BCENU_HC}, {FF_BD_HC}}};
            auto o2a = combineSharedTensors(std::move(vecL2a));
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

            // now build products blocks
            // 0) [A,F].[A,F]
            // 1) [A,F].[B,T]
            // 2) [B,F].[B,F]
            _products.push_back(TensorData{Ops::Dot::contractSingles({_data[0], false}, {_data[0], false})});
            _products.push_back(TensorData{Ops::Dot::contractSingles({_data[0], false}, {_data[1], true})});
            _products.push_back(TensorData{Ops::Dot::contractSingles({_data[1], false}, {_data[1], false})});
            // 3) [C,F].[C,F]
            // 4) [C,F].[D,T]
            // 5) [D,F].[D,F]
            _products.push_back(TensorData{Ops::Dot::contractSingles({_data[2], false}, {_data[2], false})});
            _products.push_back(TensorData{Ops::Dot::contractSingles({_data[2], false}, {_data[3], true})});
            _products.push_back(TensorData{Ops::Dot::contractSingles({_data[3], false}, {_data[3], false})});
            // 6) [A,F].[L,F]
            // 7) [B,T].[L,F]
            // 8) [L,F].[A,T]
            // 9) [L,F].[B,F]
            _products.push_back(TensorData{Ops::Dot::contractStar({_data[11], false}, {{_data[0], false}}, {{0, 0}})});
            _products.push_back(TensorData{Ops::Dot::contractStar({_data[11], false}, {{_data[1], true}}, {{0, 0}})});
            _products.push_back(TensorData{Ops::Dot::contractStar({_data[11], false}, {{_data[0], true}}, {{1, 0}})});
            _products.push_back(TensorData{Ops::Dot::contractStar({_data[11], false}, {{_data[1], false}}, {{1, 0}})});
            // 10) [A,F].[L,F].[A,T]
            // 11) [A,F].[L,F].[B,F]
            // 12) [B,T].[L,F].[A,T]
            // 13) [B,T].[L,F].[B,F]
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[11], false}, {{_data[0], false}, {_data[0], true}}, {{0, 0}, {1, 1}})});
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[11], false}, {{_data[0], false}, {_data[1], false}}, {{0, 0}, {1, 1}})});
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[11], false}, {{_data[1], true}, {_data[0], true}}, {{0, 0}, {1, 1}})});
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[11], false}, {{_data[1], true}, {_data[1], false}}, {{0, 0}, {1, 1}})});

            // 14) [A.F].[L,F].[L,T]
            // 15) [L,F].[L,F]
            auto a = TensorData{_products[6]->clone()};
            auto b = TensorData{_data[11]->clone()};
            b->conjugate();
            auto* one = static_cast<VectorContainer*>(a.get());
            auto* other = static_cast<VectorContainer*>(b.get());
            auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU_HC});
            _products.push_back(calcDot(std::move(a), *b, idxs));
            a = TensorData{_data[11]->clone()};
            b = TensorData{_data[11]->clone()};
            one = static_cast<VectorContainer*>(a.get());
            other = static_cast<VectorContainer*>(b.get());
            idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
            _products.push_back(calcDot(std::move(a), *b, idxs));

            // 16) [A, F] . [M, F] . [B, F]
            // 17) [A, T] . [M, F] . [B, T]
            // 18) [A, F] . [M, F] . [A, T]
            // 19) [C, F] . [M, F] . [D, F]
            // 20) [C, F] . [M, F] . [C, T]
            // 21) [C, F] . [A, F] . [M, F] . [B, F] . [D, F]
            // 22) [C, F] . [A, F] . [M, F] . [A, T] . [C, T]
            // 23) [A, F] . [M, F]
            // 24) [M, F] . [A, T]
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[12], false}, {{_data[0], false}, {_data[1], false}}, {{0, 0}, {2, 1}})});
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[12], false}, {{_data[1], true}, {_data[0], true}}, {{0, 0}, {2, 1}})});
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[12], false}, {{_data[0], false}, {_data[0], true}}, {{0, 0}, {2, 1}})});
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[12], false}, {{_data[2], false}, {_data[3], false}}, {{1, 0}, {3, 1}})});
            _products.push_back(TensorData{
                Ops::Dot::contractStar({_data[12], false}, {{_data[2], false}, {_data[2], true}}, {{1, 0}, {3, 1}})});
            _products.push_back(TensorData{Ops::Dot::contractStar(
                {_data[12], false}, {{_data[0], false}, {_data[2], false}, {_data[1], false}, {_data[3], false}},
                {{0, 0}, {1, 1}, {2, 2}, {3, 3}})});
            _products.push_back(TensorData{Ops::Dot::contractStar(
                {_data[12], false}, {{_data[0], false}, {_data[2], false}, {_data[0], true}, {_data[2], true}},
                {{0, 0}, {1, 1}, {2, 2}, {3, 3}})});
            _products.push_back(TensorData{Ops::Dot::contractStar({_data[12], false}, {{_data[0], false}}, {{0, 0}})});
            _products.push_back(TensorData{Ops::Dot::contractStar({_data[12], false}, {{_data[0], true}}, {{2, 0}})});

            // 25) ([C, F] . [M, F] . [D, F]) . [L, F]
            // 26) ([C, F] . [M, F] . [C, T]) . [L, F]
            a = TensorData{_products[19]->clone()};
            b = TensorData{_data[11]->clone()};
            one = static_cast<VectorContainer*>(a.get());
            other = static_cast<VectorContainer*>(b.get());
            idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
            _products.push_back(calcDot(std::move(a), *b, idxs));
            a = TensorData{_products[20]->clone()};
            b = TensorData{_data[11]->clone()};
            one = static_cast<VectorContainer*>(a.get());
            other = static_cast<VectorContainer*>(b.get());
            idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
            _products.push_back(calcDot(std::move(a), *b, idxs));

            // 27) ([A, F] . [M, F]) . [L, T] . ([A,T] . [M, F])
            a = TensorData{_products[23]->clone()};
            b = TensorData{_data[11]->clone()};
            b->conjugate();
            auto c = TensorData{_products[24]->clone()};
            one = static_cast<VectorContainer*>(a.get());
            other = static_cast<VectorContainer*>(b.get());
            idxs = one->getSameLabelPairs(*other, {WILSON_BCENU_HC});
            auto res = calcDot(std::move(a), *b, idxs);
            one = static_cast<VectorContainer*>(res.get());
            other = static_cast<VectorContainer*>(c.get());
            idxs = one->getSameLabelPairs(*other, {WILSON_BCENU});
            _products.push_back(calcDot(std::move(res), *c, idxs));
        }


        std::vector<SharedTensorData> _data;
        std::vector<TensorData> _tensors;
        std::vector<BlockIndexing> _indices;

        std::vector<TensorData> _products;
    };

    // static const IContainer* getIC(SharedTensorData& p) {
    //     return const_cast<const IContainer*>(p.get());
    // }

    static const IContainer* getIC(TensorData& p) {
        return const_cast<const IContainer*>(p.get());
    }

    // static OuterContainer* getOC(SharedTensorData& o) {
    //     return static_cast<OuterContainer*>(o.get());
    // }

    static OuterContainer* getOC(TensorData& o) {
        return static_cast<OuterContainer*>(o.get());
    }

    TEST_F(DotOuterTest, DotOneAOneB) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[1]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);


        auto val00 = _products[0]->element({});
        auto val01 = _products[1]->element({});
        EXPECT_DOUBLE_EQ(res->element({}).real(), (val01 * val00 * conj(val01 * val00)).real());
        EXPECT_DOUBLE_EQ(res->element({}).imag(), (val01 * val00 * conj(val01 * val00)).imag());
    }

    TEST_F(DotOuterTest, DotOneAOneC) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[2]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[0]->element({});
        auto val13 = conj(_products[13]->element({}));
        EXPECT_DOUBLE_EQ(res->element({}).real(), (val0 * val13 * conj(val0)).real());
        EXPECT_DOUBLE_EQ(res->element({}).imag(), (val0 * val13 * conj(val0)).imag());
    }

    TEST_F(DotOuterTest, DotOneATwoA) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[3]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[11]->element({});
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = conj(one->begin()->at(2).first->element({i3}));
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = (other->begin()->at(2).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneATwoB) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[11]->element({});
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = conj(one->begin()->at(2).first->element({i3}));
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = conj(other->begin()->at(2).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneATwoC) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[0]->element({});
        val0 *= _products[2]->element({});
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = conj(one->begin()->at(2).first->element({i3}));
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = (other->begin()->at(3).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneATwoD) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[0]->element({});
        val0 *= conj(_products[1]->element({}));
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = conj(one->begin()->at(2).first->element({i3}));
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = conj(other->begin()->at(3).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneAThree) {
        auto* one = getOC(_tensors[0]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        EXPECT_EQ(res->dims()[0], 3);
        EXPECT_EQ(res->dims()[1], 3);
        EXPECT_EQ(res->dims()[2], 3);
        EXPECT_EQ(res->dims()[3], 3);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->labels()[1], FF_BD_HC);
        EXPECT_EQ(res->labels()[2], FF_BD);
        EXPECT_EQ(res->labels()[3], FF_BD_HC);
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                auto val1 = _products[16]->element({i1, i2});
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val2 = _products[17]->element({i3, i4});
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).real(), (val1 * val2).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).imag(), (val1 * val2).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneBOneC) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[2]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[0]->element({});
        auto val1 = conj(_products[10]->element({}));
        EXPECT_DOUBLE_EQ(res->element({}).real(), (val0 * val1 * conj(val0)).real());
        EXPECT_DOUBLE_EQ(res->element({}).imag(), (val0 * val1 * conj(val0)).imag());
    }

    TEST_F(DotOuterTest, DotOneBTwoA) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[3]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[10]->element({});
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = one->begin()->at(2).first->element({i3});
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = (other->begin()->at(2).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneBTwoB) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[10]->element({});
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = one->begin()->at(2).first->element({i3});
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = conj(other->begin()->at(2).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneBTwoC) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[0]->element({});
        val0 *= conj(_products[1]->element({}));
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = one->begin()->at(2).first->element({i3});
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = (other->begin()->at(3).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneBTwoD) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[0]->element({});
        val0 *= conj(_products[0]->element({}));
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(3).first->element({i2}));
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    auto val3 = one->begin()->at(2).first->element({i3});
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val4 = conj(other->begin()->at(3).first->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posWC] = i3;
                        idx[posFFHC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneBThree) {
        auto* one = getOC(_tensors[1]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        EXPECT_EQ(res->dims()[0], 3);
        EXPECT_EQ(res->dims()[1], 3);
        EXPECT_EQ(res->dims()[2], 3);
        EXPECT_EQ(res->dims()[3], 3);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->labels()[1], FF_BD_HC);
        EXPECT_EQ(res->labels()[2], FF_BD);
        EXPECT_EQ(res->labels()[3], FF_BD_HC);
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                auto val1 = _products[18]->element({i1, i2});
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val2 = _products[18]->element({i3, i4});
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).real(), (val1 * val2).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).imag(), (val1 * val2).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneCTwoA) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[3]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(2).first->element({i2}));
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto val3 = (other->begin()->at(2).first->element({i3}));
                    for (IndexType i4 = 0; i4 < 2; ++i4) {
                        auto val4 = _products[14]->element({i4});
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posFFHC] = i3;
                        idx[posWC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneCTwoB) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(2).first->element({i2}));
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto val3 = conj(other->begin()->at(2).first->element({i3}));
                    for (IndexType i4 = 0; i4 < 2; ++i4) {
                        auto val4 = _products[14]->element({i4});
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posFFHC] = i3;
                        idx[posWC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneCTwoC) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[0]->element({});
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(2).first->element({i2}));
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto val3 = other->begin()->at(3).first->element({i3});
                    for (IndexType i4 = 0; i4 < 2; ++i4) {
                        auto val4 = conj(_products[7]->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posFFHC] = i3;
                        idx[posWC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneCTwoD) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        auto resLabels = res->labels();
        auto findLblPos = [&](IndexLabel lbl) {
            return static_cast<IndexType>(
                std::distance(resLabels.begin(), std::find(resLabels.begin(), resLabels.end(), lbl)));
        };
        IndexType posFF = findLblPos(FF_BD);
        IndexType posWCHC = findLblPos(WILSON_BCENU_HC);
        IndexType posWC = findLblPos(WILSON_BCENU);
        IndexType posFFHC = findLblPos(FF_BD_HC);
        EXPECT_LT(posFF, 4u);
        EXPECT_LT(posWCHC, 4u);
        EXPECT_LT(posWC, 4u);
        EXPECT_LT(posFFHC, 4u);
        EXPECT_EQ(res->dims()[posFF], 3);
        EXPECT_EQ(res->dims()[posWCHC], 2);
        EXPECT_EQ(res->dims()[posWC], 2);
        EXPECT_EQ(res->dims()[posFFHC], 3);
        auto val0 = _products[0]->element({});
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            auto val1 = (other->begin()->at(0).first->element({i1}));
            for (IndexType i2 = 0; i2 < 2; ++i2) {
                auto val2 = conj(one->begin()->at(2).first->element({i2}));
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    auto val3 = conj(other->begin()->at(3).first->element({i3}));
                    for (IndexType i4 = 0; i4 < 2; ++i4) {
                        auto val4 = conj(_products[6]->element({i4}));
                        IndexList idx(4);
                        idx[posFF] = i1;
                        idx[posWCHC] = i2;
                        idx[posFFHC] = i3;
                        idx[posWC] = i4;
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).real(), (val0 * val1 * val2 * val3 * val4).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element(idx).imag(), (val0 * val1 * val2 * val3 * val4).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotOneCThree) {
        auto* one = getOC(_tensors[2]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        EXPECT_EQ(res->dims()[0], 3);
        EXPECT_EQ(res->dims()[1], 3);
        EXPECT_EQ(res->dims()[2], 3);
        EXPECT_EQ(res->dims()[3], 3);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->labels()[1], FF_BD_HC);
        EXPECT_EQ(res->labels()[2], FF_BD);
        EXPECT_EQ(res->labels()[3], FF_BD_HC);
        for (IndexType i1 = 0; i1 < 3; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val0 = _products[27]->element({i3, i4, i1, i2});
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).real(), (val0).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).imag(), (val0).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotTwoATwoB) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[4]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[3]->element({});
        val0 *= _products[15]->element({});
        val0 *= conj(_products[4]->element({}));
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).real(), (val0).real());
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).imag(), (val0).imag());
    }

    TEST_F(DotOuterTest, DotTwoATwoC) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[3]->element({});
        val0 *= _products[11]->element({});
        val0 *= _products[5]->element({});
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).real(), (val0).real());
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).imag(), (val0).imag());
    }

    TEST_F(DotOuterTest, DotTwoATwoD) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[3]->element({});
        val0 *= _products[10]->element({});
        val0 *= conj(_products[4]->element({}));
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).real(), (val0).real());
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).imag(), (val0).imag());
    }

    TEST_F(DotOuterTest, DotTwoAThree) {
        auto* one = getOC(_tensors[3]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        EXPECT_EQ(res->dims()[0], 2);
        EXPECT_EQ(res->dims()[1], 3);
        EXPECT_EQ(res->dims()[2], 2);
        EXPECT_EQ(res->dims()[3], 3);
        EXPECT_EQ(res->labels()[0], WILSON_BCENU);
        EXPECT_EQ(res->labels()[1], FF_BD);
        EXPECT_EQ(res->labels()[2], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[3], FF_BD_HC);
        auto val0 = _products[25]->element({});
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val1 = _data[12]->element({i1, i2, i3, i4});
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).real(), (val0 * val1).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).imag(), (val0 * val1).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotTwoBTwoC) {
        auto* one = getOC(_tensors[4]);
        auto* other = getOC(_tensors[5]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[3]->element({});
        val0 *= _products[11]->element({});
        val0 *= conj(_products[4]->element({}));
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).real(), (val0).real());
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).imag(), (val0).imag());
    }

    TEST_F(DotOuterTest, DotTwoBTwoD) {
        auto* one = getOC(_tensors[4]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[3]->element({});
        val0 *= _products[10]->element({});
        val0 *= conj(_products[3]->element({}));
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).real(), (val0).real());
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).imag(), (val0).imag());
    }

    TEST_F(DotOuterTest, DotTwoBThree) {
        auto* one = getOC(_tensors[4]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        EXPECT_EQ(res->dims()[0], 2);
        EXPECT_EQ(res->dims()[1], 3);
        EXPECT_EQ(res->dims()[2], 2);
        EXPECT_EQ(res->dims()[3], 3);
        EXPECT_EQ(res->labels()[0], WILSON_BCENU);
        EXPECT_EQ(res->labels()[1], FF_BD);
        EXPECT_EQ(res->labels()[2], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[3], FF_BD_HC);
        auto val0 = _products[26]->element({});
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val1 = _data[12]->element({i1, i2, i3, i4});
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).real(), (val0 * val1).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).imag(), (val0 * val1).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotTwoCTwoD) {
        auto* one = getOC(_tensors[5]);
        auto* other = getOC(_tensors[6]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 0);
        auto val0 = _products[3]->element({});
        val0 *= _products[0]->element({});
        val0 *= conj(_products[1]->element({}));
        val0 *= conj(_products[4]->element({}));
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).real(), (val0).real());
        EXPECT_DOUBLE_EQ(getIC(res)->element({}).imag(), (val0).imag());
    }

    TEST_F(DotOuterTest, DotTwoCThree) {
        auto* one = getOC(_tensors[5]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        EXPECT_EQ(res->dims()[0], 2);
        EXPECT_EQ(res->dims()[1], 3);
        EXPECT_EQ(res->dims()[2], 2);
        EXPECT_EQ(res->dims()[3], 3);
        EXPECT_EQ(res->labels()[0], WILSON_BCENU);
        EXPECT_EQ(res->labels()[1], FF_BD);
        EXPECT_EQ(res->labels()[2], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[3], FF_BD_HC);
        auto val0 = _products[21]->element({});
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val1 = _data[12]->element({i1, i2, i3, i4});
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).real(), (val0 * val1).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).imag(), (val0 * val1).imag());
                    }
                }
            }
        }
    }

    TEST_F(DotOuterTest, DotTwoDThree) {
        auto* one = getOC(_tensors[6]);
        auto* other = getOC(_tensors[7]);
        auto idxs = one->getSameLabelPairs(*other, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto a = TensorData{one->clone()};
        auto b = TensorData{other->clone()};
        auto res = calcDot(std::move(a), *b, idxs);
        EXPECT_EQ(res->rank(), 4);
        EXPECT_EQ(res->dims()[0], 2);
        EXPECT_EQ(res->dims()[1], 3);
        EXPECT_EQ(res->dims()[2], 2);
        EXPECT_EQ(res->dims()[3], 3);
        EXPECT_EQ(res->labels()[0], WILSON_BCENU);
        EXPECT_EQ(res->labels()[1], FF_BD);
        EXPECT_EQ(res->labels()[2], WILSON_BCENU_HC);
        EXPECT_EQ(res->labels()[3], FF_BD_HC);
        auto val0 = _products[22]->element({});
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        auto val1 = _data[12]->element({i1, i2, i3, i4});
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).real(), (val0 * val1).real());
                        EXPECT_DOUBLE_EQ(getIC(res)->element({i1, i2, i3, i4}).imag(), (val0 * val1).imag());
                    }
                }
            }
        }
    }

} // namespace Hammer::MultiDimensional
