///
/// @file  AddAt.hh
/// @brief Sub-tensor block insertion algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_ADDAT
#define HAMMER_MATH_MULTIDIM_OPS_ADDAT

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class AddAt final {
        public:

            explicit AddAt(IndexType position, IndexType coord) noexcept;

            IContainer* operator()(VectorContainer& a, const VectorContainer& b) const;
            IContainer* operator()(SparseContainer& a, const SparseContainer& b) const;
            IContainer* operator()(VectorContainer& a, const SparseContainer& b) const;
            IContainer* operator()(SparseContainer& a, const VectorContainer& b) const;
            IContainer* operator()(SparseContainer& a, const OuterContainer& b) const;
            IContainer* operator()(VectorContainer& a, const OuterContainer& b) const;
            IContainer* operator()(OuterContainer& a, const IContainer& b) const;

            IContainer* operator()(IContainer& a, const IContainer& b) const;

            static IContainer* error(IContainer& /*unused*/, const IContainer& /*unused*/);

        private:

            IndexType _position;
            IndexType _coord;
        };

    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
