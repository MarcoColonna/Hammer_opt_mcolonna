///
/// @file  Trace.hh
/// @brief Tensor trace algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_TRACE
#define HAMMER_MATH_MULTIDIM_OPS_TRACE

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class Trace final {
        public:

            explicit Trace(IndexPairList indices) noexcept;

            IContainer* operator()(VectorContainer& a);
            IContainer* operator()(SparseContainer& a);
            IContainer* operator()(OuterContainer& a);

            static IContainer* error(IContainer& /*unused*/);

        private:

            [[nodiscard]] std::pair<IndexList, LabelsList> getNewIndexLabels(const IContainer& original) const;
            [[nodiscard]] IndexList reducedIndex(const IndexList& a) const;

            IndexPairList _indices;
            UniqueIndexList _idxSet;
        };

    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
