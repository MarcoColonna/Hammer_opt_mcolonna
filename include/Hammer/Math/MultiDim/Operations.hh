///
/// @file  Operations.hh
/// @brief Tensor operations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPERATIONS
#define HAMMER_MATH_MULTIDIM_OPERATIONS

#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    TensorData calcDot(TensorData origin, const IContainer& other, const IndexPairList& indices);
    TensorData calcTrace(TensorData origin, const IndexPairList& indices);
    TensorData calcSquare(TensorData origin);
    TensorData sum(TensorData origin, const IContainer& other);
    TensorData elementMultiply(TensorData origin, const IContainer& other);
    TensorData elementDivide(TensorData origin, const IContainer& other);
    TensorData addAt(TensorData origin, const IContainer& other, IndexType index, IndexType position);
    TensorData reOptimize(TensorData origin);
    TensorData toSparse(TensorData origin);
    TensorData toVector(TensorData origin);
    TensorData read(const Serial::FBTensor* msgreader);

    bool operator==(const IContainer& first, const IContainer& other);
    bool operator!=(const IContainer& first, const IContainer& other);
} // namespace Hammer::MultiDimensional


#endif
