///
/// @file  TestDot.cc
/// @brief Tests for Dot
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
#include "Hammer/Math/MultiDim/Ops/Dot.hh"
#include "Hammer/Math/MultiDim/Ops/Sum.hh"
#include "Hammer/Exceptions.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

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

    TEST(DotTest, SparseSparse) {
        auto t1 = makeEmptySparse({3, 2, 1}, {FF_BD, WILSON_BCENU, INTEGRATION_INDEX});
        t1->element({2, 1, 0}) = 0.1;
        t1->element({1, 0, 0}) = 0.2i;
        auto t2 = makeEmptySparse({2}, {WILSON_BCENU});
        t2->element({0}) = 0.3;
        t2->element({1}) = 0.4;
        auto res = calcDot(std::move(t1), *t2, {{1, 0}});
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 0.0);
    }

    TEST(DotTest, SparseSparse2) {
        auto t1 = makeEmptySparse({3, 2, 4, 1}, {FF_BD, WILSON_BCENU, SPIN_DSTAR, INTEGRATION_INDEX});
        t1->element({2, 0, 3, 0}) = 0.2i;
        t1->element({2, 1, 2, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.6i;
        t1->element({1, 1, 3, 0}) = 0.05;
        auto t2 = makeEmptySparse({3, 2, 7, 4}, {FF_BD_HC, WILSON_BCENU, FF_BDSTAR, SPIN_DSTAR});
        t2->element({2, 0, 5, 3}) = 0.3;
        t2->element({2, 1, 5, 2}) = 0.4;
        t2->element({2, 0, 5, 1}) = 0.8;
        t2->element({2, 1, 5, 3}) = 0.3;
        auto res = calcDot(std::move(t1), *t2, {{1, 1}, {2, 3}});
        EXPECT_EQ(res->labels().size(), 4);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->labels()[1], INTEGRATION_INDEX);
        EXPECT_EQ(res->labels()[2], FF_BD_HC);
        EXPECT_EQ(res->labels()[3], FF_BDSTAR);
        EXPECT_EQ(res->dims()[0], 3);
        EXPECT_EQ(res->dims()[1], 1);
        EXPECT_EQ(res->dims()[2], 3);
        EXPECT_EQ(res->dims()[3], 7);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0, 0}).imag(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).real(), 0.05 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).imag(), 0.6 * 0.8);
    }

    TEST(DotTest, VectorVector) {
        auto t1 = makeEmptyVector({3, 2, 1}, {FF_BD, WILSON_BCENU, INTEGRATION_INDEX});
        t1->element({2, 1, 0}) = 0.1;
        t1->element({1, 0, 0}) = 0.2i;
        auto t2 = makeEmptyVector({2}, {WILSON_BCENU});
        t2->element({0}) = 0.3;
        t2->element({1}) = 0.4;
        auto res = calcDot(std::move(t1), *t2, {{1, 0}});
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 0.0);
    }

    TEST(DotTest, VectorVector2) {
        auto t1 = makeEmptyVector({3, 2, 4, 1}, {FF_BD, WILSON_BCENU, SPIN_DSTAR, INTEGRATION_INDEX});
        t1->element({2, 0, 3, 0}) = 0.2i;
        t1->element({2, 1, 2, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.6i;
        t1->element({1, 1, 3, 0}) = 0.05;
        auto t2 = makeEmptyVector({3, 2, 7, 4}, {FF_BD_HC, WILSON_BCENU, FF_BDSTAR, SPIN_DSTAR});
        t2->element({2, 0, 5, 3}) = 0.3;
        t2->element({2, 1, 5, 2}) = 0.4;
        t2->element({2, 0, 5, 1}) = 0.8;
        t2->element({2, 1, 5, 3}) = 0.3;
        auto res = calcDot(std::move(t1), *t2, {{1, 1}, {2, 3}});
        EXPECT_EQ(res->labels().size(), 4);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->labels()[1], INTEGRATION_INDEX);
        EXPECT_EQ(res->labels()[2], FF_BD_HC);
        EXPECT_EQ(res->labels()[3], FF_BDSTAR);
        EXPECT_EQ(res->dims()[0], 3);
        EXPECT_EQ(res->dims()[1], 1);
        EXPECT_EQ(res->dims()[2], 3);
        EXPECT_EQ(res->dims()[3], 7);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({0, 0, 0, 0}).imag(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).real(), 0.05 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).imag(), 0.6 * 0.8);
    }

    TEST(DotTest, OuterOuter) {
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

        auto idxs = o1->getSameLabelPairs(*o2, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        EXPECT_EQ(idxs.size(), 2);
        EXPECT_EQ(idxs[0].first, 0);
        EXPECT_EQ(idxs[0].second, 1);
        EXPECT_EQ(idxs[1].first, 1);
        EXPECT_EQ(idxs[1].second, 3);

        auto res = calcDot(std::move(o1), *o2, idxs);

        EXPECT_EQ(res->rank(), 2);
        EXPECT_EQ(res->dims()[0], 3);
        EXPECT_EQ(res->dims()[1], 3);

        auto resEvalSingle1 =
            makeVector({3}, {FF_BD}, {1.i * (-2.) + 2. * 1., 1.i * 1.i + 2. * 5., 1.i * 3. + 2. * 1.});
        auto resEvalSingle2 =
            makeVector({3}, {FF_BD_HC}, {-2. * 3. + 1. * (-2.i), -2. * 4. + 1. * 1., -2. * (-2.) + 1. * 1.});
        vector<TensorData> vec3;
        vec3.push_back(std::move(resEvalSingle1));
        vec3.push_back(std::move(resEvalSingle2));
        auto resEvalOuter = combineTensors(vec3);
        auto resEval = makeVector({3, 3}, {FF_BD, FF_BD_HC},
                                  {-16. + 8. * 1i, -14. + 14. * 1i, 10. - 10. * 1i, -54. - 18. * 1i, -63., 45.,
                                   -6. - 22. * 1i, -14. - 21. * 1i, 10. + 15. * 1i});

        for (IndexType idx1 = 0; idx1 < 3; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).real(),
                                 getIC(resEvalOuter)->element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).imag(),
                                 getIC(resEvalOuter)->element({idx1, idx2}).imag());
                EXPECT_DOUBLE_EQ(resEval->element({idx1, idx2}).real(),
                                 getIC(resEvalOuter)->element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(resEval->element({idx1, idx2}).imag(),
                                 getIC(resEvalOuter)->element({idx1, idx2}).imag());
            }
        }
    }

    TEST(DotTest, OuterOuter2) {
        auto v1 = makeEmptyVector({6}, {WILSON_BCENU});
        v1->element({0}) = 1.i;
        v1->element({1}) = 2.;
        auto v2 = makeEmptyVector({6}, {WILSON_BCENU});
        v2->element({0}) = -2.;
        v2->element({1}) = 1.;
        auto v3 = makeEmptyVector({3, 6}, {FF_BD, WILSON_BCENU});
        v3->element({0, 0}) = -2.;
        v3->element({0, 1}) = 1.;
        v3->element({1, 0}) = 1.i;
        v3->element({1, 1}) = 5.;
        v3->element({2, 0}) = 3.;
        v3->element({2, 1}) = 1.;
        auto v4 = makeEmptyVector({3, 6}, {FF_BD, WILSON_BCENU});
        v4->element({0, 0}) = 3.;
        v4->element({0, 1}) = -2.i;
        v4->element({1, 0}) = 4.;
        v4->element({1, 1}) = 1.;
        v4->element({2, 0}) = -2.;
        v4->element({2, 1}) = 1.;

        SharedTensorData s1{v1.release()};
        SharedTensorData s2{v2.release()};
        SharedTensorData s3{v3.release()};
        SharedTensorData s4{v4.release()};

        vector<pair<SharedTensorData, bool>> vec1;
        vec1.emplace_back(s1, false);
        vec1.emplace_back(s1, true);
        vector<pair<SharedTensorData, bool>> vec2;
        vec2.emplace_back(s2, false);
        vec2.emplace_back(s2, true);
        vector<pair<SharedTensorData, bool>> vec3;
        vec3.emplace_back(s3, false);
        vec3.emplace_back(s3, true);
        vector<pair<SharedTensorData, bool>> vec4;
        vec4.emplace_back(s4, false);
        vec4.emplace_back(s4, true);

        auto o1 = combineSharedTensors(std::move(vec1));
        auto o2 = combineSharedTensors(std::move(vec2));
        auto o3 = combineSharedTensors(std::move(vec3));
        auto o4 = combineSharedTensors(std::move(vec4));

        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 0}).real(), 1.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 1}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 1}).imag(), 2.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 0}).imag(), -2.);

        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 0}).real(), 4.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 1}).real(), 1.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 1}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 1}).real(), -2.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 0}).real(), -2.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 1}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 0}).imag(), 0.);

        auto idxs = o1->getSameLabelPairs(*o3, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto idxs2 = o2->getSameLabelPairs(*o4, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});

        auto res = calcDot(std::move(o1), *o3, idxs);
        auto res2 = calcDot(std::move(o2), *o4, idxs);

        auto resEvalSingle1 = makeVector({3}, {FF_BD}, {2. * (1. - 1.i), 9, 2. + 3.i});
        SharedTensorData s5{resEvalSingle1.release()};
        auto resEvalSingle2 = makeVector({3}, {FF_BD}, {-2. * (3. + 1.i), -7., 5.});
        SharedTensorData s6{resEvalSingle2.release()};
        vector<pair<SharedTensorData, bool>> vec5;
        vector<pair<SharedTensorData, bool>> vec6;
        vec5.emplace_back(s5, false);
        vec5.emplace_back(s5, true);
        vec6.emplace_back(s6, false);
        vec6.emplace_back(s6, true);
        auto o5 = combineSharedTensors(std::move(vec5));
        auto o6 = combineSharedTensors(std::move(vec6));

        for (IndexType idx1 = 0; idx1 < 3; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).real(), getIC(o5)->element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).imag(), getIC(o5)->element({idx1, idx2}).imag());
                EXPECT_DOUBLE_EQ(getIC(res2)->element({idx1, idx2}).real(), getIC(o6)->element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(getIC(res2)->element({idx1, idx2}).imag(), getIC(o6)->element({idx1, idx2}).imag());
            }
        }
    }

    TEST(DotTest, OuterOuter3) {
        auto v1 = makeEmptySparse({6}, {WILSON_BCENU});
        v1->element({0}) = 1.i;
        v1->element({1}) = 2.;
        auto v2 = makeEmptySparse({6}, {WILSON_BCENU});
        v2->element({0}) = -2.;
        v2->element({1}) = 1.;
        auto v3 = makeEmptySparse({3, 6}, {FF_BD, WILSON_BCENU});
        v3->element({0, 0}) = -2.;
        v3->element({0, 1}) = 1.;
        v3->element({1, 0}) = 1.i;
        v3->element({1, 1}) = 5.;
        v3->element({2, 0}) = 3.;
        v3->element({2, 1}) = 1.;
        auto v4 = makeEmptySparse({3, 6}, {FF_BD, WILSON_BCENU});
        v4->element({0, 0}) = 3.;
        v4->element({0, 1}) = -2.i;
        v4->element({1, 0}) = 4.;
        v4->element({1, 1}) = 1.;
        v4->element({2, 0}) = -2.;
        v4->element({2, 1}) = 1.;

        SharedTensorData s1{v1.release()};
        SharedTensorData s2{v2.release()};
        SharedTensorData s3{v3.release()};
        SharedTensorData s4{v4.release()};

        vector<pair<SharedTensorData, bool>> vec1;
        vec1.emplace_back(s1, false);
        vec1.emplace_back(s1, true);
        vector<pair<SharedTensorData, bool>> vec2;
        vec2.emplace_back(s2, false);
        vec2.emplace_back(s2, true);
        vector<pair<SharedTensorData, bool>> vec3;
        vec3.emplace_back(s3, false);
        vec3.emplace_back(s3, true);
        vector<pair<SharedTensorData, bool>> vec4;
        vec4.emplace_back(s4, false);
        vec4.emplace_back(s4, true);

        auto o1 = combineSharedTensors(std::move(vec1));
        auto o2 = combineSharedTensors(std::move(vec2));
        auto o3 = combineSharedTensors(std::move(vec3));
        auto o4 = combineSharedTensors(std::move(vec4));

        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 0}).real(), 1.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 1}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 1}).real(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({0, 1}).imag(), 2.);
        EXPECT_DOUBLE_EQ(getOC(o1)->value({1, 0}).imag(), -2.);

        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 0}).real(), 4.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 1}).real(), 1.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 1}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 1}).real(), -2.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 0}).real(), -2.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({0, 1}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(o2)->value({1, 0}).imag(), 0.);

        auto idxs = o1->getSameLabelPairs(*o3, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto idxs2 = o2->getSameLabelPairs(*o4, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});

        auto res = calcDot(std::move(o1), *o3, idxs);
        auto res2 = calcDot(std::move(o2), *o4, idxs);

        auto resEvalSingle1 = makeVector({3}, {FF_BD}, {2. * (1. - 1.i), 9, 2. + 3.i});
        SharedTensorData s5{resEvalSingle1.release()};
        auto resEvalSingle2 = makeVector({3}, {FF_BD}, {-2. * (3. + 1.i), -7., 5.});
        SharedTensorData s6{resEvalSingle2.release()};
        vector<pair<SharedTensorData, bool>> vec5;
        vector<pair<SharedTensorData, bool>> vec6;
        vec5.emplace_back(s5, false);
        vec5.emplace_back(s5, true);
        vec6.emplace_back(s6, false);
        vec6.emplace_back(s6, true);
        auto o5 = combineSharedTensors(std::move(vec5));
        auto o6 = combineSharedTensors(std::move(vec6));

        for (IndexType idx1 = 0; idx1 < 3; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).real(), getIC(o5)->element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).imag(), getIC(o5)->element({idx1, idx2}).imag());
                EXPECT_DOUBLE_EQ(getIC(res2)->element({idx1, idx2}).real(), getIC(o6)->element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(getIC(res2)->element({idx1, idx2}).imag(), getIC(o6)->element({idx1, idx2}).imag());
            }
        }
    }

    TEST(DotTest, OuterOuter4) {
        auto v1 = makeEmptySparse({6}, {WILSON_BCENU});
        v1->element({0}) = 1.i;
        v1->element({1}) = 2.;
        auto v2 = makeEmptySparse({6}, {WILSON_BCENU});
        v2->element({0}) = -2.;
        v2->element({1}) = 1.;
        auto v3 = makeEmptySparse({3, 6}, {FF_BD, WILSON_BCENU});
        v3->element({0, 0}) = -2.;
        v3->element({0, 1}) = 1.;
        v3->element({1, 0}) = 1.i;
        v3->element({1, 1}) = 5.;
        v3->element({2, 0}) = 3.;
        v3->element({2, 1}) = 1.;
        auto v4 = makeEmptySparse({3, 6}, {FF_BD, WILSON_BCENU});
        v4->element({0, 0}) = 3.;
        v4->element({0, 1}) = -2.i;
        v4->element({1, 0}) = 4.;
        v4->element({1, 1}) = 1.;
        v4->element({2, 0}) = -2.;
        v4->element({2, 1}) = 1.;

        SharedTensorData s1{v1.release()};
        SharedTensorData s2{v2.release()};
        SharedTensorData s3{v3.release()};
        SharedTensorData s4{v4.release()};

        vector<pair<SharedTensorData, bool>> vec1;
        vec1.emplace_back(s1, false);
        vec1.emplace_back(s1, true);
        vector<pair<SharedTensorData, bool>> vec2;
        vec2.emplace_back(s2, false);
        vec2.emplace_back(s2, true);
        vector<pair<SharedTensorData, bool>> vec3;
        vec3.emplace_back(s3, false);
        vec3.emplace_back(s3, true);
        vector<pair<SharedTensorData, bool>> vec4;
        vec4.emplace_back(s4, false);
        vec4.emplace_back(s4, true);

        auto o1 = combineSharedTensors(std::move(vec1));
        auto o2 = combineSharedTensors(std::move(vec3));
        auto res1 = sum(std::move(o1), *combineSharedTensors(std::move(vec2)));
        auto res2 = sum(std::move(o2), *combineSharedTensors(std::move(vec4)));

        EXPECT_DOUBLE_EQ(getOC(res1)->value({0, 0}).real(), 5.);
        EXPECT_DOUBLE_EQ(getOC(res1)->value({1, 1}).real(), 5.);
        EXPECT_DOUBLE_EQ(getOC(res1)->value({0, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(res1)->value({1, 1}).imag(), 0.);
        EXPECT_DOUBLE_EQ(getOC(res1)->value({0, 1}).real(), -2.);
        EXPECT_DOUBLE_EQ(getOC(res1)->value({1, 0}).real(), -2.);
        EXPECT_DOUBLE_EQ(getOC(res1)->value({0, 1}).imag(), 2.);
        EXPECT_DOUBLE_EQ(getOC(res1)->value({1, 0}).imag(), -2.);

        auto idxs = res1->getSameLabelPairs(*res2, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});

        auto res = calcDot(std::move(res1), *res2, idxs);

        auto resEval = makeEmptyVector({3, 3}, {FF_BD, FF_BD_HC});
        resEval->element({0, 0}) = 74.;
        resEval->element({0, 1}) = 81. + 4.i;
        resEval->element({0, 2}) = -55. - 22.i;
        resEval->element({1, 0}) = 81. - 4.i;
        resEval->element({1, 1}) = 179.;
        resEval->element({1, 2}) = -46. - 5.i;
        resEval->element({2, 0}) = -55. + 22.i;
        resEval->element({2, 1}) = -46. + 5.i;
        resEval->element({2, 2}) = 71.;
        for (IndexType idx1 = 0; idx1 < 3; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).real(),
                                 getIC(resEval)->element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(getIC(res)->element({idx1, idx2}).imag(),
                                 getIC(resEval)->element({idx1, idx2}).imag());
            }
        }
    }

    TEST(DotTest, OuterOuter5) {
        auto v1 = makeEmptySparse({6}, {WILSON_BCENU});
        v1->element({0}) = 1.i;
        v1->element({3}) = 2.;
        auto v2 = makeEmptySparse({6}, {WILSON_BCENU});
        v2->element({0}) = -2.;
        v2->element({1}) = 1.;
        auto v3 = makeEmptySparse({3, 6}, {FF_BD, WILSON_BCENU});
        v3->element({0, 0}) = -2.;
        v3->element({0, 1}) = 1.;
        v3->element({1, 0}) = 1.i;
        v3->element({1, 1}) = 5.;
        v3->element({2, 0}) = 3.;
        v3->element({2, 1}) = 1.;
        auto v4 = makeEmptySparse({3, 6}, {FF_BD, WILSON_BCENU});
        v4->element({0, 0}) = 3.;
        v4->element({0, 1}) = -2.i;
        v4->element({1, 0}) = 4.;
        v4->element({1, 1}) = 1.;
        v4->element({2, 0}) = -2.;
        v4->element({2, 1}) = 1.;

        SharedTensorData s1{v1.release()};
        SharedTensorData s2{v2.release()};
        SharedTensorData s3{v3.release()};
        SharedTensorData s4{v4.release()};

        vector<pair<SharedTensorData, bool>> vec1;
        vec1.emplace_back(s1, false);
        vec1.emplace_back(s1, true);
        vec1.emplace_back(s1, true);
        vec1.emplace_back(s1, false);
        vector<pair<SharedTensorData, bool>> vec2;
        vec2.emplace_back(s2, false);
        vec2.emplace_back(s2, true);
        vec2.emplace_back(s2, true);
        vec2.emplace_back(s2, false);
        vector<pair<SharedTensorData, bool>> vec3;
        vec3.emplace_back(s3, false);
        vec3.emplace_back(s3, true);
        vec3.emplace_back(s3, true);
        vec3.emplace_back(s3, false);
        vector<pair<SharedTensorData, bool>> vec4;
        vec4.emplace_back(s4, false);
        vec4.emplace_back(s4, true);

        auto o1 = combineSharedTensors(std::move(vec1));
        auto res1 = sum(std::move(o1), *combineSharedTensors(std::move(vec2)));
        auto res2 = combineSharedTensors(std::move(vec3));
        auto s5 = SharedTensorData{toSparse(combineSharedTensors(std::move(vec4))).release()};
        vector<pair<SharedTensorData, bool>> vec5;
        vec5.emplace_back(s5, false);
        vec5.emplace_back(s5, true);
        auto res3 = combineSharedTensors(std::move(vec5));

        auto idxs = res1->getSameLabelPairs(*res2, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto idxs2 = res1->getSameLabelPairs(*res3, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        auto idxs3 = res2->getSameLabelPairs(*res1, {WILSON_BCENU, WILSON_BCENU_HC, FF_BD, FF_BD_HC});
        for (size_t i = 0; i < idxs.size(); ++i) {
            EXPECT_EQ(idxs[i].first, idxs2[i].first);
            EXPECT_EQ(idxs[i].second, idxs2[i].second);
        }

        auto w1 = makeVector({3}, {FF_BD}, {-2.i, -1, 3.i});
        auto w2 = makeVector({3}, {FF_BD}, {5., 5. - 2.i, -5.});
        SharedTensorData t1{w1.release()};
        SharedTensorData t2{w2.release()};
        vector<pair<SharedTensorData, bool>> wec1;
        wec1.emplace_back(t1, false);
        wec1.emplace_back(t1, true);
        wec1.emplace_back(t1, true);
        wec1.emplace_back(t1, false);
        vector<pair<SharedTensorData, bool>> wec2;
        wec2.emplace_back(t2, false);
        wec2.emplace_back(t2, true);
        wec2.emplace_back(t2, true);
        wec2.emplace_back(t2, false);
        auto tmpansA = combineSharedTensors(std::move(wec1));
        auto ansA = sum(std::move(tmpansA), *combineSharedTensors(std::move(wec2)));


        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 0, 0, 0}).real(), 17.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 3, 3, 3}).real(), 16.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 0, 0, 1}).real(), -8.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 0, 1, 0}).real(), -8.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 0, 0, 0}).real(), -8.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 1, 0, 0}).real(), -8.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 0, 1, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 0, 3, 3}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 1, 0, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 1, 1, 0}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 3, 0, 3}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 0, 0, 1}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 0, 1, 0}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 1, 0, 0}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 0, 3, 0}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 3, 0, 0}).real(), 4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 0, 0, 3}).real(), -4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 3, 3, 0}).real(), -4.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 1, 1, 1}).real(), -2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 0, 1, 1}).real(), -2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 1, 0, 1}).real(), -2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 1, 1, 0}).real(), -2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({1, 1, 1, 1}).real(), 1.);

        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 0, 0, 3}).imag(), -2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 0, 0, 0}).imag(), -2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 0, 3, 0}).imag(), 2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 3, 0, 0}).imag(), 2.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({0, 3, 3, 3}).imag(), 8.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 3, 3, 0}).imag(), 8.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 0, 3, 3}).imag(), -8.);
        EXPECT_DOUBLE_EQ(getIC(res1)->element({3, 3, 0, 3}).imag(), -8.);

        auto res1a = res1->clone();
        auto res1b = res1->clone();
        auto resA = calcDot(std::move(res1a), *res2, idxs);
        auto resB = calcDot(std::move(res1b), *res3, idxs);
        auto resC = calcDot(std::move(res2), *res1, idxs3);
        auto resD = calcDot(std::move(res3), *res1, idxs3);

        EXPECT_EQ(resA->dims().size(), 4);

        EXPECT_EQ(resB->dims().size(), 4);
        EXPECT_EQ(resC->dims().size(), 4);
        EXPECT_EQ(resD->dims().size(), 4);

        for (IndexType i1 = 0; i1 < 3; ++i1) {
            EXPECT_EQ(resA->labels()[i1], resB->labels()[i1]);
            EXPECT_EQ(resA->labels()[i1], resC->labels()[i1]);
            EXPECT_EQ(resA->labels()[i1], resD->labels()[i1]);
        }

        for (IndexType i1 = 0; i1 < 3; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 3; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        EXPECT_DOUBLE_EQ(getIC(resA)->element({i1, i2, i3, i4}).real(),
                                         getIC(ansA)->element({i1, i2, i3, i4}).real());
                        EXPECT_DOUBLE_EQ(getIC(resA)->element({i1, i2, i3, i4}).imag(),
                                         getIC(ansA)->element({i1, i2, i3, i4}).imag());
                        EXPECT_DOUBLE_EQ(getIC(resA)->element({i1, i2, i3, i4}).real(),
                                         getIC(resC)->element({i1, i2, i3, i4}).real());
                        EXPECT_DOUBLE_EQ(getIC(resA)->element({i1, i2, i3, i4}).imag(),
                                         getIC(resC)->element({i1, i2, i3, i4}).imag());
                        EXPECT_DOUBLE_EQ(getIC(resB)->element({i1, i2, i3, i4}).real(),
                                         getIC(resD)->element({i1, i2, i3, i4}).real());
                        EXPECT_DOUBLE_EQ(getIC(resB)->element({i1, i2, i3, i4}).imag(),
                                         getIC(resD)->element({i1, i2, i3, i4}).imag());
                    }
                }
            }
        }
    }


    TEST(DotTest, SparseVector) {
        auto t1 = makeEmptySparse({3, 2, 1}, {FF_BD, WILSON_BCENU, INTEGRATION_INDEX});
        t1->element({2, 1, 0}) = 0.1;
        t1->element({1, 0, 0}) = 0.2i;
        auto t2 = makeEmptyVector({2}, {WILSON_BCENU});
        t2->element({0}) = 0.3;
        t2->element({1}) = 0.4;
        auto res = calcDot(std::move(t1), *t2, {{1, 0}});
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 0.0);
    }

    TEST(DotTest, SparseVector2) {
        auto t1 = makeEmptySparse({3, 2, 4, 1}, {FF_BD, WILSON_BCENU, SPIN_DSTAR, INTEGRATION_INDEX});
        t1->element({2, 0, 3, 0}) = 0.2i;
        t1->element({2, 1, 2, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.6i;
        t1->element({1, 1, 3, 0}) = 0.05;
        auto t2 = makeEmptyVector({3, 2, 7, 4}, {FF_BD_HC, WILSON_BCENU, FF_BDSTAR, SPIN_DSTAR});
        t2->element({2, 0, 5, 3}) = 0.3;
        t2->element({2, 1, 5, 2}) = 0.4;
        t2->element({2, 0, 5, 1}) = 0.8;
        t2->element({2, 1, 5, 3}) = 0.3;
        auto res = calcDot(std::move(t1), *t2, {{1, 1}, {2, 3}});
        EXPECT_EQ(res->labels().size(), 4);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->labels()[1], INTEGRATION_INDEX);
        EXPECT_EQ(res->labels()[2], FF_BD_HC);
        EXPECT_EQ(res->labels()[3], FF_BDSTAR);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).real(), 0.05 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).imag(), 0.6 * 0.8);
    }

    TEST(DotTest, VectorSparse) {
        auto t1 = makeEmptyVector({3, 2, 1}, {FF_BD, WILSON_BCENU, INTEGRATION_INDEX});
        t1->element({2, 1, 0}) = 0.1;
        t1->element({1, 0, 0}) = 0.2i;
        auto t2 = makeEmptySparse({2}, {WILSON_BCENU});
        t2->element({0}) = 0.3;
        t2->element({1}) = 0.4;
        auto res = calcDot(std::move(t1), *t2, {{1, 0}});
        EXPECT_DOUBLE_EQ(res->element({0, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).real(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({0, 0}).imag(), 0.0);
        EXPECT_DOUBLE_EQ(res->element({1, 0}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({2, 0}).imag(), 0.0);
    }

    TEST(DotTest, VectorSparse2) {
        auto t1 = makeEmptyVector({3, 2, 4, 1}, {FF_BD, WILSON_BCENU, SPIN_DSTAR, INTEGRATION_INDEX});
        t1->element({2, 0, 3, 0}) = 0.2i;
        t1->element({2, 1, 2, 0}) = 0.1;
        t1->element({1, 0, 1, 0}) = 0.6i;
        t1->element({1, 1, 3, 0}) = 0.05;
        auto t2 = makeEmptySparse({3, 2, 7, 4}, {FF_BD_HC, WILSON_BCENU, FF_BDSTAR, SPIN_DSTAR});
        t2->element({2, 0, 5, 3}) = 0.3;
        t2->element({2, 1, 5, 2}) = 0.4;
        t2->element({2, 0, 5, 1}) = 0.8;
        t2->element({2, 1, 5, 3}) = 0.3;
        auto res = calcDot(std::move(t1), *t2, {{1, 1}, {2, 3}});
        EXPECT_EQ(res->labels().size(), 4);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->labels()[1], INTEGRATION_INDEX);
        EXPECT_EQ(res->labels()[2], FF_BD_HC);
        EXPECT_EQ(res->labels()[3], FF_BDSTAR);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).real(), 0.1 * 0.4);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 2, 5}).imag(), 0.2 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).real(), 0.05 * 0.3);
        EXPECT_DOUBLE_EQ(res->element({1, 0, 2, 5}).imag(), 0.6 * 0.8);
    }

    TEST(DotTest, OuterSparseStarTopology) {
        auto v1 = makeEmptySparse({6}, {WILSON_BCENU});
        v1->element({0}) = 1.i;
        v1->element({1}) = 2.;

        SharedTensorData s1{v1.release()};

        vector<pair<SharedTensorData, bool>> vec1;
        vec1.emplace_back(s1, false);
        vec1.emplace_back(s1, true);
        auto o1 = combineSharedTensors(std::move(vec1));

        auto t2 = makeEmptySparse({6, 6}, {WILSON_BCENU, WILSON_BCENU_HC});
        t2->element({0, 0}) = 1.;
        t2->element({1, 1}) = 1.;
        t2->element({0, 1}) = 0.5;

        auto idxs = o1->getSameLabelPairs(*t2, {WILSON_BCENU, WILSON_BCENU_HC});
        auto res = calcDot(std::move(o1), *t2, idxs);
        EXPECT_EQ(res->rank(), 0u);
    }

    TEST(DotTest, SparseOuterStarTopology) {
        auto v1 = makeEmptySparse({6}, {WILSON_BCENU});
        v1->element({0}) = 1.i;
        v1->element({1}) = 2.;
        auto v3 = makeEmptySparse({3, 6}, {FF_BD, WILSON_BCENU});
        v3->element({0, 0}) = -2.;
        v3->element({0, 1}) = 1.;
        v3->element({1, 0}) = 1.i;
        v3->element({1, 1}) = 5.;
        v3->element({2, 0}) = 3.;
        v3->element({2, 1}) = 1.;

        SharedTensorData s1{v1.release()};
        SharedTensorData s3{v3.release()};

        vector<pair<SharedTensorData, bool>> vec3;
        vec3.emplace_back(s3, false);
        vec3.emplace_back(s3, true);
        auto o3 = combineSharedTensors(std::move(vec3));

        auto t1 = makeEmptySparse({6}, {WILSON_BCENU});
        t1->element({0}) = 0.5;
        t1->element({1}) = 0.3;

        auto idxs = t1->getSameLabelPairs(*o3, {WILSON_BCENU, WILSON_BCENU_HC});
        auto res = calcDot(std::move(t1), *o3, idxs);
        EXPECT_EQ(res->rank(), 3u);
        EXPECT_EQ(res->labels()[0], FF_BD);
        EXPECT_EQ(res->dims()[0], 3);
    }

    TEST(DotTest, OuterVectorCollapseToScalar) {
        auto v1 = makeEmptyVector({2}, {WILSON_BCENU});
        v1->element({0}) = 1.i;
        v1->element({1}) = 2.;
        auto v2 = makeEmptyVector({2}, {WILSON_BCENU_HC});
        v2->element({0}) = -2.;
        v2->element({1}) = 1.;

        vector<TensorData> vec;
        vec.push_back(std::move(v1));
        vec.push_back(std::move(v2));
        auto o1 = combineTensors(std::move(vec));

        auto t2 = makeEmptyVector({2, 2}, {WILSON_BCENU, WILSON_BCENU_HC});
        t2->element({0, 0}) = 1.;
        t2->element({1, 1}) = 1.;

        auto idxs = o1->getSameLabelPairs(*t2, {WILSON_BCENU, WILSON_BCENU_HC});
        auto res = calcDot(std::move(o1), *t2, idxs);
        EXPECT_EQ(res->rank(), 0u);
    }

    TEST(DotTest, SparseSparseToScalar) {
        auto t1 = makeEmptySparse({3, 2}, {FF_BD, WILSON_BCENU});
        t1->element({1, 0}) = 0.5;
        t1->element({2, 1}) = 0.3i;
        auto t2 = makeEmptySparse({3, 2}, {FF_BD, WILSON_BCENU});
        t2->element({1, 0}) = 0.5;
        t2->element({2, 1}) = 0.3i;
        auto res = calcDot(std::move(t1), *t2, {{0, 0}, {1, 1}});
        EXPECT_EQ(res->rank(), 0u);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 0.5 * 0.5 + (0.3i * 0.3i).real());
        EXPECT_DOUBLE_EQ(res->element({}).imag(), (0.3i * 0.3i).imag());
    }

    TEST(DotTest, VectorVectorToScalar) {
        auto t1 = makeEmptyVector({3, 2}, {FF_BD, WILSON_BCENU});
        t1->element({1, 0}) = 0.5;
        t1->element({2, 1}) = 0.3i;
        auto t2 = makeEmptyVector({3, 2}, {FF_BD, WILSON_BCENU});
        t2->element({1, 0}) = 0.5;
        t2->element({2, 1}) = 0.3i;
        auto res = calcDot(std::move(t1), *t2, {{0, 0}, {1, 1}});
        EXPECT_EQ(res->rank(), 0u);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 0.5 * 0.5 + (0.3i * 0.3i).real());
        EXPECT_DOUBLE_EQ(res->element({}).imag(), (0.3i * 0.3i).imag());
    }

    TEST(DotTest, SparseVectorToScalar) {
        auto t1 = makeEmptySparse({3, 2}, {FF_BD, WILSON_BCENU});
        t1->element({1, 0}) = 0.5;
        t1->element({2, 1}) = 0.3i;
        auto t2 = makeEmptyVector({3, 2}, {FF_BD, WILSON_BCENU});
        t2->element({1, 0}) = 0.5;
        t2->element({2, 1}) = 0.3i;
        auto res = calcDot(std::move(t1), *t2, {{0, 0}, {1, 1}});
        EXPECT_EQ(res->rank(), 0u);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 0.5 * 0.5 + (0.3i * 0.3i).real());
        EXPECT_DOUBLE_EQ(res->element({}).imag(), (0.3i * 0.3i).imag());
    }

    TEST(DotTest, VectorSparseToScalar) {
        auto t1 = makeEmptyVector({3, 2}, {FF_BD, WILSON_BCENU});
        t1->element({1, 0}) = 0.5;
        t1->element({2, 1}) = 0.3i;
        auto t2 = makeEmptySparse({3, 2}, {FF_BD, WILSON_BCENU});
        t2->element({1, 0}) = 0.5;
        t2->element({2, 1}) = 0.3i;
        auto res = calcDot(std::move(t1), *t2, {{0, 0}, {1, 1}});
        EXPECT_EQ(res->rank(), 0u);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 0.5 * 0.5 + (0.3i * 0.3i).real());
        EXPECT_DOUBLE_EQ(res->element({}).imag(), (0.3i * 0.3i).imag());
    }

    TEST(DotTest, VectorOuterStarToScalar) {
        auto v1 = makeEmptyVector({3}, {WILSON_BCENU});
        v1->element({0}) = 2.;
        v1->element({1}) = 3.;
        auto v2 = makeEmptyVector({4}, {FF_BD});
        v2->element({1}) = 0.5;
        v2->element({3}) = 1.5;

        vector<TensorData> vecs;
        vecs.push_back(std::move(v1));
        vecs.push_back(std::move(v2));
        auto outerB = combineTensors(std::move(vecs));

        auto a = makeEmptyVector({3, 4}, {WILSON_BCENU, FF_BD});
        a->element({0, 1}) = 1.;
        a->element({1, 3}) = 2.;

        auto idxs = a->getSameLabelPairs(*outerB, {WILSON_BCENU, FF_BD});
        auto res = calcDot(std::move(a), *outerB, idxs);

        EXPECT_EQ(res->rank(), 0u);
        // a[0,1]*b0[0]*b1[1] + a[1,3]*b0[1]*b1[3] = 1*2*0.5 + 2*3*1.5
        EXPECT_DOUBLE_EQ(res->element({}).real(), 1. * 2. * 0.5 + 2. * 3. * 1.5);
    }

    TEST(DotTest, VectorOuterStarNonScalar) {
        auto vb1 = makeEmptyVector({3}, {WILSON_BCENU});
        vb1->element({0}) = 2.;
        vb1->element({1}) = 0.5;
        auto vb2 = makeEmptyVector({4}, {FF_BD});
        vb2->element({1}) = 1.;
        vb2->element({3}) = 3.;

        vector<TensorData> vecs;
        vecs.push_back(std::move(vb1));
        vecs.push_back(std::move(vb2));
        auto outerB = combineTensors(std::move(vecs));

        auto a = makeEmptyVector({3, 4, 5}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        a->element({0, 1, 2}) = 1.;
        a->element({1, 3, 4}) = 2.i;

        auto idxs = a->getSameLabelPairs(*outerB, {WILSON_BCENU, FF_BD});
        auto res = calcDot(std::move(a), *outerB, idxs);

        EXPECT_EQ(res->rank(), 1u);
        EXPECT_EQ(res->dims()[0], 5u);
        // res[2] = a[0,1,2]*b0[0]*b1[1] = 1*2*1 = 2
        EXPECT_DOUBLE_EQ(res->element({2}).real(), 2.);
        EXPECT_DOUBLE_EQ(res->element({2}).imag(), 0.);
        // res[4] = a[1,3,4]*b0[1]*b1[3] = 2i*0.5*3 = 3i
        EXPECT_DOUBLE_EQ(res->element({4}).real(), 0.);
        EXPECT_DOUBLE_EQ(res->element({4}).imag(), 3.);
    }

    TEST(DotTest, VectorOuterBoomerang) {
        // Outer with two rank-2 sub-tensors; Vector contracted with one index from each.
        // Result has same rank as a but dims at contracted positions are replaced by free dims of b sub-tensors.
        auto vb0 = makeEmptyVector({3, 5}, {WILSON_BCENU, SPIN_DSTAR});
        vb0->element({0, 2}) = 1.;
        vb0->element({1, 3}) = 2.;
        auto vb1 = makeEmptyVector({4, 6}, {FF_BD, FF_BDSTAR});
        vb1->element({1, 0}) = 0.5;
        vb1->element({3, 4}) = 1.5;

        vector<TensorData> vecs;
        vecs.push_back(std::move(vb0));
        vecs.push_back(std::move(vb1));
        auto outerB = combineTensors(std::move(vecs));

        auto a = makeEmptyVector({3, 4, 2}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        a->element({0, 1, 0}) = 3.;
        a->element({1, 3, 1}) = 4.;

        auto idxs = a->getSameLabelPairs(*outerB, {WILSON_BCENU, FF_BD, SPIN_DSTAR, FF_BDSTAR});
        EXPECT_EQ(idxs.size(), 2u);
        auto res = calcDot(std::move(a), *outerB, idxs);

        EXPECT_EQ(res->rank(), 3u);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 0}).real(), 1.5);
        EXPECT_DOUBLE_EQ(res->element({3, 4, 1}).real(), 12.);
    }

    TEST(DotTest, SparseOuterBoomerang) {
        auto vb0 = makeEmptySparse({3, 5}, {WILSON_BCENU, SPIN_DSTAR});
        vb0->element({0, 2}) = 1.;
        vb0->element({1, 3}) = 2.;
        auto vb1 = makeEmptySparse({4, 6}, {FF_BD, FF_BDSTAR});
        vb1->element({1, 0}) = 0.5;
        vb1->element({3, 4}) = 1.5;

        vector<TensorData> vecs;
        vecs.push_back(std::move(vb0));
        vecs.push_back(std::move(vb1));
        auto outerB = combineTensors(std::move(vecs));

        auto a = makeEmptySparse({3, 4, 2}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX});
        a->element({0, 1, 0}) = 3.;
        a->element({1, 3, 1}) = 4.;

        auto idxs = a->getSameLabelPairs(*outerB, {WILSON_BCENU, FF_BD, SPIN_DSTAR, FF_BDSTAR});
        EXPECT_EQ(idxs.size(), 2u);
        auto res = calcDot(std::move(a), *outerB, idxs);

        EXPECT_EQ(res->rank(), 3u);
        EXPECT_DOUBLE_EQ(res->element({2, 0, 0}).real(), 1.5);
        EXPECT_DOUBLE_EQ(res->element({3, 4, 1}).real(), 12.);
    }

    TEST(DotTest, OuterVectorPartialCollapse) {
        auto va = makeEmptyVector({3}, {WILSON_BCENU});
        va->element({0}) = 2.;
        va->element({2}) = -1.;
        auto vb = makeEmptyVector({4}, {FF_BD});
        vb->element({1}) = 0.5;
        auto vc = makeEmptyVector({5}, {INTEGRATION_INDEX});
        vc->element({3}) = 1.5;

        vector<TensorData> vecs;
        vecs.push_back(std::move(va));
        vecs.push_back(std::move(vb));
        vecs.push_back(std::move(vc));
        auto outerA = combineTensors(std::move(vecs));

        auto b = makeEmptyVector({3}, {WILSON_BCENU});
        b->element({0}) = 1.;
        b->element({2}) = 3.;

        auto idxs = outerA->getSameLabelPairs(*b, {WILSON_BCENU});
        EXPECT_EQ(idxs.size(), 1u);
        auto res = calcDot(std::move(outerA), *b, idxs);

        EXPECT_EQ(res->rank(), 2u);
        const IContainer* cres = res.get();
        EXPECT_DOUBLE_EQ(cres->element({1, 3}).real(), -0.75);
    }

    TEST(DotTest, OuterVectorNonStarGeneral) {
        auto vb0 = makeEmptyVector({3, 4}, {WILSON_BCENU, FF_BD});
        vb0->element({0, 1}) = 2.;
        vb0->element({1, 2}) = 3.;
        auto vb1 = makeEmptyVector({5, 6}, {INTEGRATION_INDEX, SPIN_DSTAR});
        vb1->element({2, 3}) = 0.5;
        vb1->element({3, 4}) = 1.;

        vector<TensorData> vecs;
        vecs.push_back(std::move(vb0));
        vecs.push_back(std::move(vb1));
        auto outerA = combineTensors(std::move(vecs));

        auto b = makeEmptyVector({3, 4, 5, 6}, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX, SPIN_DSTAR});
        b->element({0, 1, 2, 3}) = 1.;
        b->element({1, 2, 3, 4}) = 2.;

        auto idxs = outerA->getSameLabelPairs(*b, {WILSON_BCENU, FF_BD, INTEGRATION_INDEX, SPIN_DSTAR});
        EXPECT_EQ(idxs.size(), 4u);
        auto res = calcDot(std::move(outerA), *b, idxs);

        EXPECT_EQ(res->rank(), 0u);
        EXPECT_DOUBLE_EQ(res->element({}).real(), 7.);
    }

    TEST(DotTest, ErrorPath) {
        auto scal = makeEmptyScalar();
        scal->element({}) = 1.;
        auto t2 = makeEmptyVector({2}, {WILSON_BCENU});
        t2->element({0}) = 1.;
        EXPECT_THROW(calcDot(std::move(scal), *t2, {}), Hammer::Error);
    }

} // namespace Hammer::MultiDimensional
