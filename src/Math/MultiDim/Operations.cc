///
/// @file  Operations.cc
/// @brief Tensor operations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/Operations.hh"
#include "Hammer/Math/MultiDim/OperationDefs.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Math/MultiDim/Ops/Dot.hh"
#include "Hammer/Math/MultiDim/Ops/Trace.hh"
#include "Hammer/Math/MultiDim/Ops/Sum.hh"
#include "Hammer/Math/MultiDim/Ops/AddAt.hh"
#include "Hammer/Math/MultiDim/Ops/Multiply.hh"
#include "Hammer/Math/MultiDim/Ops/Divide.hh"
#include "Hammer/Math/MultiDim/Ops/OuterSquare.hh"
#include "Hammer/Math/MultiDim/Ops/Optimize.hh"
#include "Hammer/Math/MultiDim/Ops/Convert.hh"
#include "Hammer/Math/MultiDim/Ops/Compare.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"


using namespace std;

namespace Hammer::MultiDimensional {

    TensorData calcDot(TensorData origin, const IContainer& other, const IndexPairList& indices) {
        Ops::Dot dotter{indices};
        return calc2(std::move(origin), other, dotter, "dot");
    }

    TensorData calcTrace(TensorData origin, const IndexPairList& indices) {
        Ops::Trace tracer{indices};
        return calc1(std::move(origin), tracer, "trace");
    }

    TensorData calcSquare(TensorData origin) {
        Ops::OuterSquare squarer{};
        return calc1(std::move(origin), squarer, "square");
    }

    TensorData sum(TensorData origin, const IContainer& other) {
        Ops::Sum summer{};
        return calc2(std::move(origin), other, summer, "sum");
    }

    TensorData elementMultiply(TensorData origin, const IContainer& other) {
        Ops::Multiply multiplier{};
        return calc2(std::move(origin), other, multiplier, "element-multiply");
    }

    TensorData elementDivide(TensorData origin, const IContainer& other) {
        Ops::Divide divider{};
        return calc2(std::move(origin), other, divider, "element-divide");
    }

    TensorData addAt(TensorData origin, const IContainer& other, IndexType index, IndexType position) {
        Ops::AddAt inserter{index, position};
        return calc2(std::move(origin), other, inserter, "add-at");
    }

    TensorData read(const Serial::FBTensor* msgreader) {
        TensorData result;
        switch (msgreader->data_type()) {
        case Serial::FBTensorTypes::FBComplex: {
            result.reset(static_cast<IContainer*>(new ScalarContainer{msgreader->data_as_FBComplex()}));
            break;
        }
        case Serial::FBTensorTypes::FBSingleTensor: {
            const auto* data = msgreader->data_as_FBSingleTensor();
            if (data->sparse()) {
                result.reset(static_cast<IContainer*>(new SparseContainer{data}));
            } else {
                result.reset(static_cast<IContainer*>(new VectorContainer{data}));
            }
            break;
        }
        case Serial::FBTensorTypes::FBTensorList: {
            const auto* data = msgreader->data_as_FBTensorList();
            result.reset(static_cast<IContainer*>(new OuterContainer{data}));
            break;
        }
        case Serial::FBTensorTypes::NONE:
            break;
        }
        return result;
    }

    TensorData reOptimize(TensorData origin) {
        Ops::Optimize optimizer{};
        return calc1(std::move(origin), optimizer, "optimize");
    }

    TensorData toSparse(TensorData origin) {
        Ops::Convert converter{true};
        return calc1(std::move(origin), converter, "convert to sparse");
    }

    TensorData toVector(TensorData origin) {
        Ops::Convert converter{false};
        return calc1(std::move(origin), converter, "convert to vector");
    }

    bool operator==(const IContainer& first, const IContainer& other) {
        Ops::Compare comparer{};
        bool result = ComparisonDispatch<Ops::Compare>::Go(first, other, comparer);
        return result;
    }

    bool operator!=(const IContainer& first, const IContainer& other) {
        return !(first == other);
    }


} // namespace Hammer::MultiDimensional
