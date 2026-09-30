///
/// @file  OperationDefs.hh
/// @brief Tensor operations helper functions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <string>

#include "Hammer/Tools/Loki.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Math/MultiDim/Operations.fhh"

namespace Hammer::MultiDimensional {

    class IContainer;

    template <class Executor>
    using BinaryDispatch = Loki::StaticDoubleDispatcher<Executor, IContainer, TensorDataTypes, false, const IContainer,
                                                        ConstTensorDataTypes, IContainer*>;
    template <class Executor>
    using UnaryDispatch = Loki::StaticSingleDispatcher<Executor, IContainer, TensorDataTypes, IContainer*>;

    template <class Executor>
    using ComparisonDispatch = Loki::StaticDoubleDispatcher<Executor, const IContainer, ConstTensorDataTypes, false,
                                                            const IContainer, ConstTensorDataTypes, bool>;

    template <class Ops, class TensorPtr>
    TensorPtr calc2(TensorPtr origin, const IContainer& other, Ops op, const std::string& opName) {
        IContainer* result = BinaryDispatch<Ops>::Go(*origin, other, op);
        if (result == nullptr) {
            throw Error("Unable to " + opName);
        }
        if (result != origin.get()) {
            origin.reset(result);
        }
        return origin;
    }

    template <class Ops, class TensorPtr>
    TensorPtr calc1(TensorPtr origin, Ops op, const std::string& opName) {
        IContainer* result = UnaryDispatch<Ops>::Go(*origin, op);
        if (result == nullptr) {
            throw Error("Unable to " + opName);
        }
        if (result != origin.get()) {
            origin.reset(result);
        }
        return origin;
    }
} // namespace Hammer::MultiDimensional
