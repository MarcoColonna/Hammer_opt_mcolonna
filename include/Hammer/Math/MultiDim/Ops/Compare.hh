///
/// @file  Compare.hh
/// @brief Tensor storage type conversion algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_Compare
#define HAMMER_MATH_MULTIDIM_OPS_Compare

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class Compare final {
        public:

            Compare() noexcept = default;

            static bool error(const IContainer& /*unused*/, const IContainer& /*unused*/);

            bool operator()(const VectorContainer& a, const VectorContainer& b);
            bool operator()(const SparseContainer& a, const SparseContainer& b);
            bool operator()(const VectorContainer& a, const SparseContainer& b);
            bool operator()(const SparseContainer& a, const VectorContainer& b);

            bool operator()(const IContainer& a, const IContainer& b);

        private:

            static bool checkStructure(const IContainer& a, const IContainer& b);
        };


    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
