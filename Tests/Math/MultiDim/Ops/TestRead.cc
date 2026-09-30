///
/// @file  TestRead.cc
/// @brief Tests for Read
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
#include "Hammer/Exceptions.hh"
#include "Hammer/Tools/HammerSerial.hh"

#include "gtest/gtest.h"

using namespace std;

namespace Hammer::MultiDimensional {

    TEST(ReadTest, Sparse) {
        auto t1 = makeEmptySparse({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({1, 2, 0}) = 0.1;
        t1->element({0, 1, 0}) = 0.2 + 0.3i;
        auto builder = make_unique<flatbuffers::FlatBufferBuilder>(16ul * 1024ul);
        auto res = t1->write(builder.get());
        Serial::FBTensorBuilder serialTensor{*builder};
        serialTensor.add_data_type(res.second);
        serialTensor.add_data(res.first);
        auto sert = serialTensor.Finish();
        builder->Finish(sert);
        auto* start = builder->GetBufferPointer();
        const auto* rtensor = flatbuffers::GetRoot<Serial::FBTensor>(start);
        auto res2 = read(rtensor);
        EXPECT_TRUE(t1->compare(*res2));
    }

    TEST(ReadTest, Vector) {
        auto t1 = makeEmptyVector({2, 3, 1}, {WILSON_BCENU, FF_BD, SPIN_DSTAR});
        t1->element({1, 0, 0}) = 0.1;
        t1->element({1, 1, 0}) = 0.2;
        t1->element({1, 2, 0}) = 0.3;
        t1->element({0, 0, 0}) = 0.2 + 0.3i;
        t1->element({0, 1, 0}) = 0.4 + 0.3i;
        t1->element({0, 2, 0}) = 0.6 + 0.3i;
        auto builder = make_unique<flatbuffers::FlatBufferBuilder>(16ul * 1024ul);
        auto res = t1->write(builder.get());
        Serial::FBTensorBuilder serialTensor{*builder};
        serialTensor.add_data_type(res.second);
        serialTensor.add_data(res.first);
        auto sert = serialTensor.Finish();
        builder->Finish(sert);
        auto* start = builder->GetBufferPointer();
        const auto* rtensor = flatbuffers::GetRoot<Serial::FBTensor>(start);
        auto res2 = read(rtensor);
        EXPECT_TRUE(t1->compare(*res2));
    }

    TEST(ReadTest, Outer) {
        auto t = makeEmptyVector({2}, {WILSON_BCENU});
        t->element({1}) = 0.1;
        t->element({0}) = 0.2i;
        auto t2 = makeOuterSquare(t);
        auto builder = make_unique<flatbuffers::FlatBufferBuilder>(16ul * 1024ul);
        auto res = t2->write(builder.get());
        Serial::FBTensorBuilder serialTensor{*builder};
        serialTensor.add_data_type(res.second);
        serialTensor.add_data(res.first);
        auto sert = serialTensor.Finish();
        builder->Finish(sert);
        auto* start = builder->GetBufferPointer();
        const auto* rtensor = flatbuffers::GetRoot<Serial::FBTensor>(start);
        auto res2 = read(rtensor);
        EXPECT_TRUE(t2->compare(*res2));
    }

    TEST(ReadTest, Outer2) {
        auto t = makeEmptyVector({2}, {WILSON_BCENU});
        t->element({1}) = 0.1;
        t->element({0}) = 0.2i;
        auto t2 = makeEmptyVector({2}, {WILSON_BCENU});
        t2->element({1}) = 0.7;
        t2->element({0}) = 0.9i;
        auto t3 = makeOuterSquare(t);
        auto t4 = makeOuterSquare(t2);
        auto tmp = sum(std::move(t3), *t4);
        auto builder = make_unique<flatbuffers::FlatBufferBuilder>(16ul * 1024ul);
        auto res = tmp->write(builder.get());
        Serial::FBTensorBuilder serialTensor{*builder};
        serialTensor.add_data_type(res.second);
        serialTensor.add_data(res.first);
        auto sert = serialTensor.Finish();
        builder->Finish(sert);
        auto* start = builder->GetBufferPointer();
        const auto* rtensor = flatbuffers::GetRoot<Serial::FBTensor>(start);
        auto res2 = read(rtensor);
        EXPECT_TRUE(tmp->compare(*res2));
    }

} // namespace Hammer::MultiDimensional
