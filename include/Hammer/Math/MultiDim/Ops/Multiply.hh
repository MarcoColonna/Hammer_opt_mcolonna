///
/// @file  Multiply.hh
/// @brief Tensor element-wise multiplication algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_MULTIPLY
#define HAMMER_MATH_MULTIDIM_OPS_MULTIPLY

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class Multiply final {
        public:

            Multiply() noexcept = default;

            IContainer* operator()(VectorContainer& a, const VectorContainer& b);
            IContainer* operator()(SparseContainer& a, const SparseContainer& b);
            IContainer* operator()(VectorContainer& a, const SparseContainer& b);
            IContainer* operator()(SparseContainer& a, const VectorContainer& b);
            IContainer* operator()(OuterContainer& a, const VectorContainer& b);
            IContainer* operator()(OuterContainer& a, const SparseContainer& b);
            IContainer* operator()(OuterContainer& a, const OuterContainer& b);

            IContainer* operator()(IContainer& a, const IContainer& b);

            static IContainer* error(IContainer& /*unused*/, const IContainer& /*unused*/);
        };

    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
