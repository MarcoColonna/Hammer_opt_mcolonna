///
/// @file  TestOuterContainer.cc
/// @brief Tests for OuterContainer
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    TEST(OuterTest, ConstructionAccess) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2;
        const auto t2 = makeOuterSquare(t);
        EXPECT_EQ(t2->rank(), 6);
        EXPECT_EQ(t->rank(), 3);
        EXPECT_EQ(t2->dims()[3], 2);
        EXPECT_EQ(t2->dims()[4], 3);
        EXPECT_EQ(t2->dims()[5], 1);
        EXPECT_EQ(t2->labels()[3], WILSON_BCENU_HC);
        EXPECT_EQ(t2->labels()[4], FF_BD_HC);
        EXPECT_EQ(t2->labels()[5], SPIN_DSTAR_HC);
        t->element({0, 0, 0}) = 0.5i;
        const auto t3 = makeOuterSquare(std::move(t));
        EXPECT_EQ(t3->rank(), 6);
        EXPECT_EQ(t3->dims()[3], 2);
        EXPECT_EQ(t3->dims()[4], 3);
        EXPECT_EQ(t3->dims()[5], 1);
        EXPECT_EQ(t3->labels()[3], WILSON_BCENU_HC);
        EXPECT_EQ(t3->labels()[4], FF_BD_HC);
        EXPECT_EQ(t3->labels()[5], SPIN_DSTAR_HC);
        EXPECT_FALSE(t);
        /// @todo fix this const crap
        const OuterContainer* o2 = static_cast<OuterContainer*>(t2.get());
        const OuterContainer* o3 = static_cast<OuterContainer*>(t3.get());
        EXPECT_DOUBLE_EQ(o2->element({0, 0, 0, 0, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(o3->element({0, 0, 0, 0, 0, 0}).real(), 0.5 * 0.5);
        EXPECT_DOUBLE_EQ(o3->element({0, 0, 0, 0, 0, 0}).imag(), 0.);
    }

    TEST(OuterTest, ConstructionAccess2) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2;
        auto t2 = makeEmptySparse({5, 4, 2}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t2->element({3, 3, 0}) = 0.1;
        t2->element({2, 1, 1}) = 0.2i;
        vector<TensorData> vec;
        vec.push_back(std::move(t));
        vec.push_back(std::move(t2));
        auto res = combineTensors(vec);
        EXPECT_EQ(res->rank(), 6);
        EXPECT_EQ(vec[0]->rank(), 3);
        EXPECT_EQ(vec[1]->rank(), 3);
        EXPECT_EQ(res->dims()[3], 5);
        EXPECT_EQ(res->dims()[4], 4);
        EXPECT_EQ(res->dims()[5], 2);
        EXPECT_EQ(res->labels()[3], WILSON_BCENU);
        EXPECT_EQ(res->labels()[4], FF_BD);
        EXPECT_EQ(res->labels()[5], SPIN_DSTAR);
        vec[0]->element({0, 0, 0}) = 0.5;
        auto res2 = combineTensors(std::move(vec));
        EXPECT_EQ(res2->rank(), 6);
        EXPECT_FALSE(vec[0]); // NOLINT(bugprone-use-after-move) -- intentional: verifies ownership was transferred
        EXPECT_FALSE(vec[1]); // NOLINT(bugprone-use-after-move)
        EXPECT_EQ(res2->dims()[3], 5);
        EXPECT_EQ(res2->dims()[4], 4);
        EXPECT_EQ(res2->dims()[5], 2);
        EXPECT_EQ(res2->labels()[3], WILSON_BCENU);
        EXPECT_EQ(res2->labels()[4], FF_BD);
        EXPECT_EQ(res2->labels()[5], SPIN_DSTAR);
        /// @todo fix this const crap
        const OuterContainer* o2 = static_cast<OuterContainer*>(res.get());
        const OuterContainer* o3 = static_cast<OuterContainer*>(res2.get());
        EXPECT_DOUBLE_EQ(o2->element({0, 0, 0, 0, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(o3->element({0, 0, 0, 0, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(o2->element({0, 1, 0, 0, 0, 0}).real(), 0.);
        EXPECT_DOUBLE_EQ(o3->element({0, 0, 0, 2, 1, 1}).imag(), 0.5 * 0.2);
    }

    TEST(OuterTest, ValueAccess) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.3 + 0.2i;
        auto t2 = makeOuterSquare(t);
        unique_ptr<OuterContainer> o2{static_cast<OuterContainer*>(t2.release())};
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 1, 2, 0}).real(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 1, 2, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 0, 1, 0}).real(), 0.3 * 0.3 + 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 0, 1, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 0, 1, 0}).real(), 0.1 * 0.3);
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 0, 1, 0}).imag(), -0.1 * 0.2);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 1, 2, 0}).real(), 0.3 * 0.1);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 1, 2, 0}).imag(), 0.2 * 0.1);
    }

    TEST(OuterTest, CopyClone) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1 + 0.2i;
        t->element({0, 1, 0}) = 0.2;
        auto t2 = makeOuterSquare(t);
        auto t4 = t2->clone();
        EXPECT_EQ(t2->rank(), t4->rank());
        EXPECT_EQ(t2->dims()[0], t4->dims()[0]);
        EXPECT_EQ(t2->dims()[1], t4->dims()[1]);
        EXPECT_EQ(t2->dims()[2], t4->dims()[2]);
        EXPECT_EQ(t2->dims()[3], t4->dims()[3]);
        EXPECT_EQ(t2->dims()[4], t4->dims()[4]);
        EXPECT_EQ(t2->dims()[5], t4->dims()[5]);
        EXPECT_EQ(t2->labels()[0], t4->labels()[0]);
        EXPECT_EQ(t2->labels()[1], t4->labels()[1]);
        EXPECT_EQ(t2->labels()[2], t4->labels()[2]);
        EXPECT_EQ(t2->labels()[3], t4->labels()[3]);
        EXPECT_EQ(t2->labels()[4], t4->labels()[4]);
        EXPECT_EQ(t2->labels()[5], t4->labels()[5]);
        /// @todo fix this const crap
        const OuterContainer* o2 = static_cast<OuterContainer*>(t2.get());
        const OuterContainer* o4 = static_cast<OuterContainer*>(t4.get());
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        EXPECT_DOUBLE_EQ(o2->element({i1, i2, 0, i3, i4, 0}).real(),
                                         o4->element({i1, i2, 0, i3, i4, 0}).real());
                        EXPECT_DOUBLE_EQ(o2->element({i1, i2, 0, i3, i4, 0}).imag(),
                                         o4->element({i1, i2, 0, i3, i4, 0}).imag());
                    }
                }
            }
        }
        t->element({0, 0, 0}) = 0.5;
        auto t3 = makeOuterSquare(std::move(t));
        unique_ptr<const OuterContainer> o3{static_cast<OuterContainer*>(t3.release())};
        const auto& t5 = *o3;
        EXPECT_EQ(t5.rank(), t2->rank());
        EXPECT_EQ(t5.dims()[0], o3->dims()[0]);
        EXPECT_EQ(t5.dims()[1], o3->dims()[1]);
        EXPECT_EQ(t5.dims()[2], o3->dims()[2]);
        EXPECT_EQ(t5.labels()[0], o3->labels()[0]);
        EXPECT_EQ(t5.labels()[1], o3->labels()[1]);
        EXPECT_EQ(t5.labels()[2], o3->labels()[2]);
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        EXPECT_DOUBLE_EQ(t5.element({i1, i2, 0, i3, i4, 0}).real(),
                                         o3->element({i1, i2, 0, i3, i4, 0}).real());
                        EXPECT_DOUBLE_EQ(t5.element({i1, i2, 0, i3, i4, 0}).imag(),
                                         o3->element({i1, i2, 0, i3, i4, 0}).imag());
                    }
                }
            }
        }
    }

    TEST(OuterTest, Comparison) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        auto t2 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        auto t3 = makeEmptySparse({4, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        auto t4 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD_HC, SPIN_DSTAR});
        auto t5 = makeEmptyScalar();
        auto o1 = makeOuterSquare(t1);
        auto o2 = makeOuterSquare(t2);
        auto o3 = makeOuterSquare(t3);
        auto o4 = makeOuterSquare(t4);
        auto o5 = makeOuterSquare(t5);
        EXPECT_TRUE(o1->compare(*o2));
        EXPECT_FALSE(o1->compare(*o3));
        EXPECT_FALSE(o1->compare(*o4));
        EXPECT_FALSE(o1->compare(*o5));
        /// @todo fix this const crap
        const OuterContainer* res1 = static_cast<OuterContainer*>(o1.get());
        res1->element({0, 1, 0, 0, 1, 0}) = 0.3;
        EXPECT_FALSE(o1->compare(*o2));
    }

    TEST(OuterTest, ComponentMultiply) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2;
        auto res = makeOuterSquare(t);
        /// @todo fix this const crap
        const OuterContainer* o = static_cast<OuterContainer*>(res.get());
        EXPECT_DOUBLE_EQ(o->element({1, 2, 0, 1, 2, 0}).real(), 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(o->element({0, 1, 0, 0, 1, 0}).real(), 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(o->element({0, 0, 0, 0, 0, 0}).real(), 0.0);
        (*res) *= 3.;
        EXPECT_DOUBLE_EQ(o->element({1, 2, 0, 1, 2, 0}).real(), 3 * 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(o->element({0, 1, 0, 0, 1, 0}).real(), 3 * 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(o->element({0, 0, 0, 0, 0, 0}).real(), 0.0);
    }

    TEST(OuterTest, Conjugate) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1 + 0.2i;
        t->element({0, 1, 0}) = 0.2;
        auto t2 = makeOuterSquare(t);
        t2->conjugate();
        unique_ptr<OuterContainer> o2{static_cast<OuterContainer*>(t2.release())};
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 1, 2, 0}).real(), 0.1 * 0.1 + 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 1, 2, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 0, 1, 0}).real(), 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 0, 1, 0}).imag(), 0.);
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 0, 1, 0}).real(), 0.1 * 0.2);
        EXPECT_DOUBLE_EQ(o2->value({1, 2, 0, 0, 1, 0}).imag(), -0.2 * 0.2);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 1, 2, 0}).real(), 0.2 * 0.1);
        EXPECT_DOUBLE_EQ(o2->value({0, 1, 0, 1, 2, 0}).imag(), 0.2 * 0.2);
        EXPECT_EQ(o2->labels()[0], WILSON_BCENU_HC);
        EXPECT_EQ(o2->labels()[1], FF_BD_HC);
        EXPECT_EQ(o2->labels()[2], SPIN_DSTAR_HC);
        EXPECT_EQ(o2->labels()[3], WILSON_BCENU);
        EXPECT_EQ(o2->labels()[4], FF_BD);
        EXPECT_EQ(o2->labels()[5], SPIN_DSTAR);
    }


    TEST(OuterTest, IsSameShape) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        auto t2 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        auto t3 = makeEmptySparse({4, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        auto t4 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD_HC, SPIN_DSTAR});
        auto o1 = makeOuterSquare(t1);
        auto o2 = makeOuterSquare(t2);
        auto o3 = makeOuterSquare(t3);
        auto o4 = makeOuterSquare(t4);
        EXPECT_TRUE(o1->isSameShape(*o2));
        EXPECT_FALSE(o1->isSameShape(*o3));
        EXPECT_FALSE(o1->isSameShape(*o4));
        EXPECT_FALSE(o1->canAddAt(*t1, WILSON_BCENU, 0));
    }

    TEST(OuterTest, NumAddends) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        auto t2 = makeEmptySparse({5, 4, 2}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t2->element({3, 3, 0}) = 0.2i;
        std::vector<TensorData> vec;
        vec.push_back(std::move(t));
        vec.push_back(std::move(t2));
        auto res = combineTensors(std::move(vec));
        const OuterContainer* oc = static_cast<OuterContainer*>(res.get());
        EXPECT_EQ(oc->numAddends(), 1);
        EXPECT_EQ(oc->entrySize(), sizeof(std::complex<double>));
        EXPECT_GT(oc->dataSize(), 0ul);
    }

    TEST(OuterTest, ValueIterator) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1 + 0.2i;
        t->element({0, 1, 0}) = 0.3;
        auto t2 = makeOuterSquare(t);
        const OuterContainer* oc = static_cast<OuterContainer*>(t2.get());
        IndexList idx{1, 2, 0, 1, 2, 0};
        auto v1 = oc->value(idx.begin(), idx.end());
        EXPECT_DOUBLE_EQ(v1.real(), 0.1 * 0.1 + 0.2 * 0.2);
        EXPECT_DOUBLE_EQ(v1.imag(), 0.);
        std::vector<IndexList> splitIdx{{1, 2, 0}, {1, 2, 0}};
        auto v2 = oc->value(splitIdx);
        EXPECT_DOUBLE_EQ(v2.real(), v1.real());
        EXPECT_DOUBLE_EQ(v2.imag(), v1.imag());
    }

    TEST(OuterTest, BeginEndIterators) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.5;
        auto t2 = makeEmptySparse({5, 4, 2}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t2->element({2, 1, 1}) = 0.3i;
        std::vector<TensorData> vec;
        vec.push_back(std::move(t));
        vec.push_back(std::move(t2));
        auto res = combineTensors(std::move(vec));
        auto* oc = static_cast<OuterContainer*>(res.get());
        size_t count = 0;
        for (auto it = oc->begin(); it != oc->end(); ++it) {
            ++count;
        }
        EXPECT_EQ(count, 1ul);
        const OuterContainer* coc = oc;
        count = 0;
        for (auto it = coc->begin(); it != coc->end(); ++it) {
            ++count;
        }
        EXPECT_EQ(count, 1ul);
    }

    TEST(OuterTest, ComplexMultiply) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        t->element({0, 1, 0}) = 0.2;
        auto res = makeOuterSquare(t);
        const OuterContainer* o = static_cast<OuterContainer*>(res.get());
        EXPECT_DOUBLE_EQ(o->element({1, 2, 0, 1, 2, 0}).real(), 0.1 * 0.1);
        std::complex<double> factor{2., 0.};
        (*res) *= factor;
        EXPECT_DOUBLE_EQ(o->element({1, 2, 0, 1, 2, 0}).real(), 2. * 0.1 * 0.1);
        EXPECT_DOUBLE_EQ(o->element({0, 1, 0, 0, 1, 0}).real(), 2. * 0.2 * 0.2);
    }

    TEST(OuterTest, ComplexMultiplyImagFactor) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({0, 1, 0}) = 0.5;
        auto res = makeOuterSquare(t);
        const OuterContainer* o = static_cast<OuterContainer*>(res.get());
        EXPECT_DOUBLE_EQ(o->element({0, 1, 0, 0, 1, 0}).real(), 0.25);
        std::complex<double> factor{0., 1.};
        (*res) *= factor;
        // each sub-tensor gets multiplied by sqrt(i) = (1+i)/sqrt(2);
        // result element = [0.5*(1+i)/sqrt(2)]^2 = 0.25 * (2i/2) = 0.25i
        EXPECT_NEAR(o->element({0, 1, 0, 0, 1, 0}).real(), 0., 1e-12);
        EXPECT_NEAR(o->element({0, 1, 0, 0, 1, 0}).imag(), 0.25, 1e-12);
    }

    TEST(OuterTest, SharedDataMultiply) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        SharedTensorData s{t.release()};
        std::vector<std::pair<SharedTensorData, bool>> vec;
        vec.emplace_back(s, false);
        vec.emplace_back(s, true);
        auto res = combineSharedTensors(std::move(vec));
        const OuterContainer* o = static_cast<OuterContainer*>(res.get());
        EXPECT_DOUBLE_EQ(o->element({1, 2, 0, 1, 2, 0}).real(), 0.1 * 0.1);
        (*res) *= 4.;
        EXPECT_DOUBLE_EQ(o->element({1, 2, 0, 1, 2, 0}).real(), 4. * 0.1 * 0.1);
    }

    TEST(OuterTest, CombineEmpty) {
        std::vector<TensorData> empty;
        auto res = combineTensors(empty);
        EXPECT_EQ(res->rank(), 0ul);
        std::vector<TensorData> empty2;
        auto res2 = combineTensors(std::move(empty2));
        EXPECT_EQ(res2->rank(), 0ul);
        std::vector<std::pair<SharedTensorData, bool>> empty3;
        auto res3 = combineSharedTensors(std::move(empty3));
        EXPECT_EQ(res3->rank(), 0ul);
    }

    TEST(OuterTest, CompareDifferentSharedDataFlag) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({0, 1, 0}) = 0.3;
        auto o1 = makeOuterSquare(t1);
        SharedTensorData s{t1.release()};
        std::vector<std::pair<SharedTensorData, bool>> vec;
        vec.emplace_back(s, false);
        vec.emplace_back(s, true);
        auto o2 = combineSharedTensors(std::move(vec));
        EXPECT_FALSE(o1->compare(*o2));
    }

    TEST(OuterTest, HasNaNs) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        auto t2 = makeOuterSquare(t);
        const OuterContainer* oc = static_cast<OuterContainer*>(t2.get());
        EXPECT_FALSE(oc->hasNaNs());
    }

    TEST(OuterTest, SharedDataConstructorConjugated) {
        auto t = makeEmptySparse({2, 3}, {WILSON_BCENU, FF_BD});
        t->element({1, 2}) = 0.3 + 0.4i;
        SharedTensorData s{t.release()};
        OuterContainer oc{s, true};
        EXPECT_EQ(oc.rank(), 4u);
        EXPECT_EQ(oc.labels()[0], WILSON_BCENU);
        EXPECT_EQ(oc.labels()[1], FF_BD);
        EXPECT_EQ(oc.labels()[2], WILSON_BCENU_HC);
        EXPECT_EQ(oc.labels()[3], FF_BD_HC);
        const OuterContainer& coc = oc;
        EXPECT_DOUBLE_EQ(coc.element({1, 2, 1, 2}).real(), 0.3 * 0.3 + 0.4 * 0.4);
        EXPECT_NEAR(coc.element({1, 2, 1, 2}).imag(), 0., 1e-14);
    }

    TEST(OuterTest, AssignmentOperatorSharedData) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        SharedTensorData s{t.release()};
        std::vector<std::pair<SharedTensorData, bool>> vec;
        vec.emplace_back(s, false);
        vec.emplace_back(s, true);
        auto shared_res = combineSharedTensors(std::move(vec));
        auto* o1 = static_cast<OuterContainer*>(shared_res.get());
        auto nonshared_base = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        auto nonshared_res = makeOuterSquare(nonshared_base);
        auto* o2 = static_cast<OuterContainer*>(nonshared_res.get());
        *o2 = *o1;
        EXPECT_EQ(o2->rank(), o1->rank());
        const OuterContainer* co1 = o1;
        const OuterContainer* co2 = o2;
        EXPECT_DOUBLE_EQ(co2->element({1, 2, 0, 1, 2, 0}).real(), co1->element({1, 2, 0, 1, 2, 0}).real());
    }

    TEST(OuterTest, CopyConstructorSharedData) {
        auto t = makeEmptySparse({2, 3}, {WILSON_BCENU, FF_BD});
        t->element({1, 2}) = 0.5;
        SharedTensorData s{t.release()};
        std::vector<std::pair<SharedTensorData, bool>> vec;
        vec.emplace_back(s, false);
        vec.emplace_back(s, true);
        auto res = combineSharedTensors(std::move(vec));
        const OuterContainer* src = static_cast<OuterContainer*>(res.get());
        OuterContainer copy{*src};
        EXPECT_EQ(copy.rank(), src->rank());
        const OuterContainer& ccopy = copy;
        EXPECT_DOUBLE_EQ(ccopy.element({1, 2, 1, 2}).real(), src->element({1, 2, 1, 2}).real());
    }

    TEST(OuterTest, ClearMethod) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1;
        auto res = makeOuterSquare(t);
        auto* oc = static_cast<OuterContainer*>(res.get());
        EXPECT_EQ(oc->numAddends(), 1u);
        oc->clear();
        EXPECT_EQ(oc->numAddends(), 0u);
    }

    TEST(OuterTest, AssignmentOperator) {
        auto t = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t->element({1, 2, 0}) = 0.1 + 0.2i;
        t->element({0, 1, 0}) = 0.3;
        auto t2 = makeOuterSquare(t);
        auto t3 = makeOuterSquare(t);
        unique_ptr<OuterContainer> o2{static_cast<OuterContainer*>(t2.release())};
        unique_ptr<OuterContainer> o3{static_cast<OuterContainer*>(t3.release())};
        *o3 = *o2;
        EXPECT_EQ(o3->rank(), o2->rank());
        const OuterContainer* co2 = o2.get();
        const OuterContainer* co3 = o3.get();
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        EXPECT_DOUBLE_EQ(co3->element({i1, i2, 0, i3, i4, 0}).real(),
                                         co2->element({i1, i2, 0, i3, i4, 0}).real());
                        EXPECT_DOUBLE_EQ(co3->element({i1, i2, 0, i3, i4, 0}).imag(),
                                         co2->element({i1, i2, 0, i3, i4, 0}).imag());
                    }
                }
            }
        }
    }

} // namespace Hammer::MultiDimensional
