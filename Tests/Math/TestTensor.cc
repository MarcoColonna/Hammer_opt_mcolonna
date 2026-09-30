///
/// @file  TestTensor.cc
/// @brief Tests for Tensor
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Math/Tensor.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Tools/HammerSerial.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    TEST(TensorTest, Construction) {
        Tensor t{"Trial", MD::makeEmptyVector({2, 3, 1}, {SPIN_TAUP, SPIN_MUP, NONE})};
        EXPECT_EQ(t.rank(), 3);
        EXPECT_EQ(t.name(), "Trial");
        EXPECT_EQ(t.dims()[0], 2);
        EXPECT_EQ(t.dims()[1], 3);
        EXPECT_EQ(t.dims()[2], 1);
        EXPECT_EQ(t.labels()[0], SPIN_TAUP);
        EXPECT_EQ(t.labels()[1], SPIN_MUP);
        EXPECT_EQ(t.labels()[2], NONE);
    }

    TEST(TensorTest, SetGetValues) {
        Tensor t{"Trial", MD::makeVector({2, 3, 1}, {SPIN_TAUP, SPIN_MUP, NONE}, {1., 0., 0., -1., 0., 1.i})};
        EXPECT_NEAR(t.value(0, 0, 0).real(), 1.0, 1.e-5);
        EXPECT_NEAR(t.value(0, 0, 0).imag(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(0, 1, 0).real(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(0, 1, 0).imag(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(0, 2, 0).real(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(0, 2, 0).imag(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 0, 0).real(), -1.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 0, 0).imag(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 1, 0).real(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 1, 0).imag(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 2, 0).real(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 2, 0).imag(), 1.0, 1.e-5);
        t.setValue(3. + 2.i, 1, 2, 0);
        EXPECT_NEAR(t.value(1, 2, 0).real(), 3.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 2, 0).imag(), 2.0, 1.e-5);
    }

    TEST(TensorTest, Clear) {
        Tensor t{"Trial", MD::makeVector({2, 3, 1}, {SPIN_TAUP, SPIN_MUP, NONE}, {1., 0., 0., -1., 0., 1.i})};
        EXPECT_NEAR(t.value(1, 0, 0).real(), -1.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 2, 0).imag(), 1.0, 1.e-5);
        EXPECT_NEAR(t.value(0, 0, 0).real(), 1.0, 1.e-5);
        t.clearData();
        EXPECT_NEAR(t.value(1, 0, 0).real(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(1, 2, 0).imag(), 0.0, 1.e-5);
        EXPECT_NEAR(t.value(0, 0, 0).real(), 0.0, 1.e-5);
        EXPECT_EQ(t.rank(), 3);
        EXPECT_EQ(t.dims()[0], 2);
        EXPECT_EQ(t.dims()[1], 3);
        EXPECT_EQ(t.dims()[2], 1);
        EXPECT_EQ(t.labels()[0], SPIN_TAUP);
        EXPECT_EQ(t.labels()[1], SPIN_MUP);
        EXPECT_EQ(t.labels()[2], NONE);
        EXPECT_EQ(t.name(), "Trial");
    }

    TEST(TensorTest, Trace) {
        Tensor t{"Trial1", MD::makeVector({2, 3, 2}, {SPIN_TAUP, SPIN_MUP, SPIN_TAUP_HC},
                                          {1., 0., 0., -1., 0., 1.i, -1.i, 0., 0., 1., 1., 0.})};
        t.spinSum();
        EXPECT_EQ(t.rank(), 1);
        EXPECT_EQ(t.dims()[0], 3);
        EXPECT_EQ(t.labels()[0], SPIN_MUP);
        EXPECT_NEAR(0., t.element({0}).imag(), 1.e-5);
        EXPECT_NEAR(1., t.element({0}).real(), 1.e-5);
        EXPECT_NEAR(0., t.element({1}).imag(), 1.e-5);
        EXPECT_NEAR(1., t.element({1}).real(), 1.e-5);
        EXPECT_NEAR(0., t.element({2}).imag(), 1.e-5);
        EXPECT_NEAR(0., t.element({2}).real(), 1.e-5);
    }

    TEST(TensorTest, RankZero) {
        Tensor t{"Trial1", MD::makeVector({2, 2}, {SPIN_TAUP, SPIN_TAUP_HC}, {1., 0., 0., -1.})};
        t.spinSum();
        EXPECT_EQ(t.rank(), 0);
        EXPECT_EQ(t.labels().size(), 0);
        EXPECT_NEAR(0., t.element().imag(), 1.e-5);
        EXPECT_NEAR(0., t.element().real(), 1.e-5);
        // Tensor t2;
        // t2.setDimensions({2, 3, 2});
        // t2.setLabels({SPIN_TAUP, SPIN_TAUP_HC});
        // t2.setName("Trial2");
    }

    TEST(TensorTest, spinAverage) {
        Tensor t{"Trial1", MD::makeVector({2, 1, 2, 2, 1, 2},
                                          {SPIN_TAUP, SPIN_NUTAU_REF, WILSON_BCTAUNU, SPIN_TAUP_HC, SPIN_NUTAU_HC_REF,
                                           WILSON_BCTAUNU_HC},
                                          {1., 0., 0., -1., 0., 1.i, -1.i, 0., 0., 1., 1., 0., 2., 5., 4., 3.})};
        EXPECT_EQ(t.rank(), 6);
        t.spinAverage();
        EXPECT_EQ(t.rank(), 2);
        EXPECT_EQ(t.dims()[0], 2);
        EXPECT_EQ(t.labels()[0], WILSON_BCTAUNU);
        EXPECT_NEAR(0., t.element({0, 0}).imag(), 1.e-5);
        EXPECT_NEAR(2. / 4., t.element({0, 0}).real(), 1.e-5);
        EXPECT_NEAR(0., t.element({0, 1}).imag(), 1.e-5);
        EXPECT_NEAR(0. / 4., t.element({0, 1}).real(), 1.e-5);
        EXPECT_NEAR(0., t.element({1, 0}).imag(), 1.e-5);
        EXPECT_NEAR(4. / 4., t.element({1, 0}).real(), 1.e-5);
        EXPECT_NEAR(1. / 4., t.element({1, 1}).imag(), 1.e-5);
        EXPECT_NEAR(3. / 4., t.element({1, 1}).real(), 1.e-5);
    }

    TEST(TensorTest, ScalarProduct) {
        Tensor t{"Trial1", MD::makeVector({2, 3, 1}, {SPIN_TAUP, SPIN_MUP, NONE}, {1., 0., 0., -1., 0., 1.i})};
        Tensor t2 = 2. * t;
        Tensor t3 = t * 3.5;
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                EXPECT_NEAR(2. * t.element({idx1, idx2, 0}).real(), t2.element({idx1, idx2, 0}).real(), 1.e-5);
                EXPECT_NEAR(2. * t.element({idx1, idx2, 0}).imag(), t2.element({idx1, idx2, 0}).imag(), 1.e-5);
                EXPECT_NEAR(3.5 * t.element({idx1, idx2, 0}).real(), t3.element({idx1, idx2, 0}).real(), 1.e-5);
                EXPECT_NEAR(3.5 * t.element({idx1, idx2, 0}).imag(), t3.element({idx1, idx2, 0}).imag(), 1.e-5);
            }
        }
    }

    TEST(TensorTest, DotProduct) {
        Tensor t{"Trial1", MD::makeVector({2, 2}, {SPIN_TAUP, SPIN_MUP}, {1., 1., 2., -1.})};
        Tensor t2{"Trial2", MD::makeVector({2, 2}, {SPIN_TAUP, SPIN_GAMMA}, {0., 1.i, 2., -1.i})};
        Tensor t3{"Trial3", MD::makeVector({2, 2}, {SPIN_MUP, SPIN_TAUP}, {1.i, 1., 5.i, -1.})};
        Tensor r1 = dot(t, t2);
        Tensor r2 = dot(t, t3, {SPIN_MUP});
        Tensor r3 = dot(t, t3, {SPIN_MUP, SPIN_TAUP});
        EXPECT_EQ(r1.rank(), 2);
        EXPECT_EQ(r2.rank(), 2);
        EXPECT_EQ(r3.rank(), 0);
        EXPECT_EQ(r1.dims()[0], 2);
        EXPECT_EQ(r2.dims()[0], 2);
        EXPECT_EQ(r1.dims()[1], 2);
        EXPECT_EQ(r2.dims()[1], 2);
        EXPECT_EQ(r1.labels().size(), 2);
        EXPECT_EQ(r1.labels()[0], SPIN_MUP);
        EXPECT_EQ(r1.labels()[1], SPIN_GAMMA);
        EXPECT_EQ(r2.labels().size(), 2);
        EXPECT_EQ(r2.labels()[0], SPIN_TAUP);
        EXPECT_EQ(r2.labels()[1], SPIN_TAUP);
        EXPECT_NEAR(r1.element({0, 0}).real(), 4., 1.e-5);
        EXPECT_NEAR(r1.element({0, 0}).imag(), 0., 1.e-5);
        EXPECT_NEAR(r1.element({0, 1}).real(), 0., 1.e-5);
        EXPECT_NEAR(r1.element({0, 1}).imag(), -1., 1.e-5);
        EXPECT_NEAR(r1.element({1, 0}).real(), -2., 1.e-5);
        EXPECT_NEAR(r1.element({1, 0}).imag(), 0., 1.e-5);
        EXPECT_NEAR(r1.element({1, 1}).real(), 0., 1.e-5);
        EXPECT_NEAR(r1.element({1, 1}).imag(), 2., 1.e-5);
        EXPECT_NEAR(r2.element({0, 0}).real(), 0., 1.e-5);
        EXPECT_NEAR(r2.element({0, 0}).imag(), 6., 1.e-5);
        EXPECT_NEAR(r2.element({0, 1}).real(), 0., 1.e-5);
        EXPECT_NEAR(r2.element({0, 1}).imag(), 0., 1.e-5);
        EXPECT_NEAR(r2.element({1, 0}).real(), 0., 1.e-5);
        EXPECT_NEAR(r2.element({1, 0}).imag(), -3., 1.e-5);
        EXPECT_NEAR(r2.element({1, 1}).real(), 3., 1.e-5);
        EXPECT_NEAR(r2.element({1, 1}).imag(), 0., 1.e-5);
        EXPECT_NEAR(r3.element().real(), 3., 1.e-5);
        EXPECT_NEAR(r3.element().imag(), 6., 1.e-5);
    }

    TEST(TensorTest, DotProduct_2) {
        Tensor t{"Trial1",
                 MD::makeVector({2, 6}, {SPIN_TAUP, FF_BD}, {2., 2., 2., 2., 2., 2., 1., 1., 1., 1., 1., 1.})};
        Tensor t2 = t;
        Tensor t1{"Trial1", MD::makeVector({6}, {FF_BD}, {1., 2., 3., 4., 5., 6.})};

        t.dot(t1);

        EXPECT_NEAR(t.element({0}).real(), 42., 1e-5);
        EXPECT_NEAR(t.element({1}).real(), 21., 1e-5);

        t1.dot(t2);

        EXPECT_NEAR(t1.element({0}).real(), 42., 1e-5);
        EXPECT_NEAR(t1.element({1}).real(), 21., 1e-5);
    }


    TEST(TensorTest, DotProduct_3) {
        Tensor t{"Trial1", MD::makeEmptyVector({2, 2, 2, 11, 8},
                                               {SPIN_NUTAU, SPIN_NUTAU_REF, SPIN_TAUP, WILSON_BCTAUNU, FF_BDSTAR})};

        t.element({0, 0, 0, 0, 1}) = 1.0;
        t.element({0, 0, 0, 0, 2}) = 2.0;
        t.element({0, 0, 0, 0, 3}) = 3.0;
        t.element({0, 0, 0, 0, 4}) = 4.0;
        t.element({0, 0, 1, 0, 1}) = 5.0;
        t.element({0, 0, 1, 0, 2}) = 6.0;
        t.element({0, 0, 1, 0, 4}) = 7.0;

        EXPECT_FALSE(isZero(t.element({0, 0, 0, 0, 1})));
        EXPECT_TRUE(isZero(t.element({0, 0, 0, 0, 0})));

        Tensor t1{"Trial1", MD::makeEmptyVector({2, 2}, {SPIN_NUTAU_REF, SPIN_TAUP})};

        t1.element({0, 0}) = 0.3;
        t1.element({1, 1}) = 0.7;
        t1.element({1, 0}) = 1.1;
        t1.element({0, 1}) = 1.3;


        t.dot(t1);

        EXPECT_NEAR(t.element({0, 0, 0}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({0, 0, 1}).real(), 1. * 0.3 + 5.0 * 1.3, 1e-5);
        EXPECT_NEAR(t.element({0, 0, 2}).real(), 2. * 0.3 + 6.0 * 1.3, 1e-5);
        EXPECT_NEAR(t.element({0, 0, 3}).real(), 3. * 0.3, 1e-5);
        EXPECT_NEAR(t.element({0, 0, 4}).real(), 4. * 0.3 + 7.0 * 1.3, 1e-5);
        EXPECT_NEAR(t.element({0, 0, 5}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({0, 0, 6}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({0, 0, 7}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 0}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 1}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 2}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 3}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 4}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 5}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 6}).real(), 0., 1e-5);
        EXPECT_NEAR(t.element({1, 0, 7}).real(), 0., 1e-5);
    }

    TEST(TensorTest, OuterProduct) {
        Tensor t{"Trial1", MD::makeVector({2, 3}, {SPIN_TAUP, SPIN_MUP}, {1., 3.5i, -2.5i, -7., -3., 15.})};
        const Tensor r = outerSquare(t);
        EXPECT_EQ(r.rank(), 4);
        EXPECT_EQ(r.dims()[0], 2);
        EXPECT_EQ(r.dims()[1], 3);
        EXPECT_EQ(r.dims()[2], 2);
        EXPECT_EQ(r.dims()[3], 3);
        EXPECT_EQ(r.labels()[0], SPIN_TAUP);
        EXPECT_EQ(r.labels()[1], SPIN_MUP);
        EXPECT_EQ(r.labels()[2], SPIN_TAUP_HC);
        EXPECT_EQ(r.labels()[3], SPIN_MUP_HC);
        for (IndexType i1 = 0; i1 < 2; ++i1) {
            for (IndexType i2 = 0; i2 < 3; ++i2) {
                for (IndexType i3 = 0; i3 < 2; ++i3) {
                    for (IndexType i4 = 0; i4 < 3; ++i4) {
                        EXPECT_NEAR(r.element({i1, i2, i3, i4}).real(),
                                    (t.value(i1, i2) * std::conj(t.value(i3, i4))).real(), 1.e-3);
                        EXPECT_NEAR(r.element({i1, i2, i3, i4}).imag(),
                                    (t.value(i1, i2) * std::conj(t.value(i3, i4))).imag(), 1.e-3);
                    }
                }
            }
        }
    }

    TEST(TensorTest, AddAt) {
        Tensor t{"AddAt", MD::makeEmptyVector({2, 3, 2}, {FF_BD, FF_BD_HC, INTEGRATION_INDEX})};

        Tensor t1{"AddAt1", MD::makeVector({2, 3}, {FF_BD, FF_BD_HC}, {1., 0., 0., -1., 0., 1.i})};
        Tensor t2{"AddAt2", MD::makeVector({2, 3}, {FF_BD, FF_BD_HC}, {5., 6., 6., 1. + 6. * 1i, 6. * 1i, 5. * 1i})};

        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                for (IndexType idx3 = 0; idx3 < 2; ++idx3) {
                    EXPECT_DOUBLE_EQ(t.element({idx1, idx2, idx3}).real(), 0.);
                    EXPECT_DOUBLE_EQ(t.element({idx1, idx2, idx3}).imag(), 0.);
                }
            }
        }
        t.addAt(t1, INTEGRATION_INDEX, 0);
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 0}).real(), t1.element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 0}).imag(), t1.element({idx1, idx2}).imag());
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 1}).real(), 0.);
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 1}).imag(), 0.);
            }
        }
        t.addAt(t2, INTEGRATION_INDEX, 1);
        for (IndexType idx1 = 0; idx1 < 2; ++idx1) {
            for (IndexType idx2 = 0; idx2 < 3; ++idx2) {
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 0}).real(), t1.element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 0}).imag(), t1.element({idx1, idx2}).imag());
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 1}).real(), t2.element({idx1, idx2}).real());
                EXPECT_DOUBLE_EQ(t.element({idx1, idx2, 1}).imag(), t2.element({idx1, idx2}).imag());
            }
        }
    }

    TEST(TensorTest, Serialization) {
        Tensor t1{"AddAt1", MD::makeVector({2, 3}, {FF_BD, FF_BD_HC}, {1., 0., 0., -1., 0., 1.i})};
        Tensor t2{"AddAt2", MD::makeVector({2, 3}, {FF_BD, FF_BD_HC}, {5., 6., 6., 1. + 6. * 1i, 6. * 1i, 5. * 1i})};

        auto builder = make_unique<flatbuffers::FlatBufferBuilder>();
        flatbuffers::Offset<Serial::FBTensor> msg;
        builder->Clear();
        t1.write(builder.get(), &msg);
        builder->Finish(msg);
        uint8_t* buf = builder->GetBufferPointer();
        // unsigned int size = builder->GetSize();
        const auto* ser1 = flatbuffers::GetRoot<Serial::FBTensor>(buf);
        EXPECT_NE(ser1, nullptr);
        if (ser1 != nullptr) {
            Tensor t3;
            t3.read(ser1);
            EXPECT_EQ(t3.rank(), t1.rank());
            EXPECT_NE(t3.name(), t1.name());
            EXPECT_EQ(t3.name(), "Tensor");
            for (size_t idx = 0; idx < min(t1.rank(), t3.rank()); ++idx) {
                EXPECT_EQ(t1.dims()[idx], t3.dims()[idx]);
                EXPECT_EQ(t1.labels()[idx], t3.labels()[idx]);
            }
            for (IndexType idx1 = 0; idx1 < min(t1.dims()[0], t3.dims()[0]); ++idx1) {
                for (IndexType idx2 = 0; idx2 < min(t1.dims()[1], t3.dims()[1]); ++idx2) {
                    EXPECT_EQ(t1.element({idx1, idx2}).real(), t3.element({idx1, idx2}).real());
                    EXPECT_EQ(t1.element({idx1, idx2}).imag(), t3.element({idx1, idx2}).imag());
                }
            }
        }
        builder->Clear();
        t2.write(builder.get(), &msg);
        builder->Finish(msg);
        buf = builder->GetBufferPointer();
        // size = builder->GetSize();
        const auto* ser2 = flatbuffers::GetRoot<Serial::FBTensor>(buf);
        EXPECT_NE(ser2, nullptr);
        if (ser2 != nullptr) {
            Tensor t4;
            t4.read(ser2);
            EXPECT_EQ(t4.rank(), t2.rank());
            EXPECT_NE(t4.name(), t2.name());
            EXPECT_EQ(t4.name(), "Tensor");
            for (size_t idx = 0; idx < min(t2.rank(), t4.rank()); ++idx) {
                EXPECT_EQ(t2.dims()[idx], t4.dims()[idx]);
                EXPECT_EQ(t2.labels()[idx], t4.labels()[idx]);
            }
            for (IndexType idx1 = 0; idx1 < min(t2.dims()[0], t4.dims()[0]); ++idx1) {
                for (IndexType idx2 = 0; idx2 < min(t2.dims()[1], t4.dims()[1]); ++idx2) {
                    EXPECT_EQ(t2.element({idx1, idx2}).real(), t4.element({idx1, idx2}).real());
                    EXPECT_EQ(t2.element({idx1, idx2}).imag(), t4.element({idx1, idx2}).imag());
                }
            }
        }
    }

    TEST(TensorTest, HasLabels) {
        Tensor twc{"WC", MD::makeVector({2}, {WILSON_BCTAUNU}, {1., 2.})};
        EXPECT_TRUE(twc.hasWCLabels());
        EXPECT_FALSE(twc.hasFFLabels());
        EXPECT_FALSE(twc.hasFFVarLabels());

        Tensor tff{"FF", MD::makeVector({6}, {FF_BD}, {1., 2., 3., 4., 5., 6.})};
        EXPECT_FALSE(tff.hasWCLabels());
        EXPECT_TRUE(tff.hasFFLabels());
        EXPECT_FALSE(tff.hasFFVarLabels());

        Tensor tffv{"FFVar", MD::makeVector({2}, {FF_BD_VAR}, {1., 2.})};
        EXPECT_FALSE(tffv.hasWCLabels());
        EXPECT_FALSE(tffv.hasFFLabels());
        EXPECT_TRUE(tffv.hasFFVarLabels());

        Tensor tspin{"Spin", MD::makeVector({2}, {SPIN_TAUP}, {1., 2.})};
        EXPECT_FALSE(tspin.hasWCLabels());
        EXPECT_FALSE(tspin.hasFFLabels());
        EXPECT_FALSE(tspin.hasFFVarLabels());
    }

    TEST(TensorTest, InPlaceOperators) {
        Tensor t{"T", MD::makeVector({2}, {SPIN_TAUP}, {2., 3.})};
        t *= 2.0;
        EXPECT_NEAR(t.value(0).real(), 4., 1e-10);
        EXPECT_NEAR(t.value(1).real(), 6., 1e-10);
        EXPECT_NEAR(t.value(0).imag(), 0., 1e-10);

        Tensor t2{"T2", MD::makeVector({2}, {SPIN_TAUP}, {1., 2.})};
        t2 *= std::complex<double>{0., 1.};
        EXPECT_NEAR(t2.value(0).real(), 0., 1e-10);
        EXPECT_NEAR(t2.value(0).imag(), 1., 1e-10);
        EXPECT_NEAR(t2.value(1).real(), 0., 1e-10);
        EXPECT_NEAR(t2.value(1).imag(), 2., 1e-10);

        Tensor ta{"TA", MD::makeVector({2}, {SPIN_TAUP}, {1., 2.})};
        Tensor tb{"TB", MD::makeVector({2}, {SPIN_TAUP}, {3., 4.})};
        ta += tb;
        EXPECT_NEAR(ta.value(0).real(), 4., 1e-10);
        EXPECT_NEAR(ta.value(1).real(), 6., 1e-10);
    }

    TEST(TensorTest, ElementwiseOps) {
        Tensor ta{"TA", MD::makeVector({3}, {SPIN_TAUP}, {1., 2., 3.})};
        Tensor tb{"TB", MD::makeVector({3}, {SPIN_TAUP}, {2., 3., 4.})};
        ta.elementMultiplyBy(tb);
        EXPECT_NEAR(ta.value(0).real(), 2., 1e-10);
        EXPECT_NEAR(ta.value(1).real(), 6., 1e-10);
        EXPECT_NEAR(ta.value(2).real(), 12., 1e-10);

        Tensor tc{"TC", MD::makeVector({3}, {SPIN_TAUP}, {2., 6., 12.})};
        Tensor td{"TD", MD::makeVector({3}, {SPIN_TAUP}, {2., 3., 4.})};
        tc.elementDivideBy(td);
        EXPECT_NEAR(tc.value(0).real(), 1., 1e-10);
        EXPECT_NEAR(tc.value(1).real(), 2., 1e-10);
        EXPECT_NEAR(tc.value(2).real(), 3., 1e-10);
    }

    TEST(TensorTest, IsEqualToAndShape) {
        Tensor ta{"TA", MD::makeVector({2}, {SPIN_TAUP}, {1., 2.})};
        Tensor tb{"TB", MD::makeVector({2}, {SPIN_TAUP}, {1., 2.})};
        Tensor tc{"TC", MD::makeVector({2}, {SPIN_TAUP}, {1., 3.})};

        EXPECT_TRUE(ta.isEqualTo(tb));
        EXPECT_FALSE(ta.isEqualTo(tc));

        EXPECT_TRUE(ta.isSameLabelShape(tb));
        EXPECT_TRUE(ta.isSameLabelShape(tc));
    }

    TEST(TensorTest, HasNaNs) {
        Tensor t1{"T1", MD::makeVector({2}, {SPIN_TAUP}, {1., 2.})};
        EXPECT_FALSE(t1.hasNaNs());

        Tensor t2{"T2", MD::makeVector({2}, {SPIN_TAUP}, {std::numeric_limits<double>::quiet_NaN(), 2.})};
        EXPECT_TRUE(t2.hasNaNs());
    }

    TEST(TensorTest, ToVectorSparse) {
        Tensor t{"T", MD::makeVector({2, 3}, {SPIN_TAUP, SPIN_MUP}, {1., 0., 2., 0., 3., 0.})};
        t.toVector();
        EXPECT_EQ(t.rank(), 2);
        EXPECT_NEAR(t.element({0, 0}).real(), 1., 1e-10);
        EXPECT_NEAR(t.element({0, 2}).real(), 2., 1e-10);
        EXPECT_NEAR(t.element({1, 1}).real(), 3., 1e-10);

        Tensor ts{"TS", MD::makeVector({2, 3}, {SPIN_TAUP, SPIN_MUP}, {1., 0., 2., 0., 3., 0.})};
        ts.toSparse();
        EXPECT_EQ(ts.rank(), 2);
        EXPECT_NEAR(ts.element({0, 0}).real(), 1., 1e-10);
        EXPECT_NEAR(ts.element({0, 2}).real(), 2., 1e-10);
        EXPECT_NEAR(ts.element({1, 1}).real(), 3., 1e-10);
        EXPECT_NEAR(ts.element({0, 1}).real(), 0., 1e-10); // zero entry preserved
    }

} // namespace Hammer
