///
/// @file  Optimize.hh
/// @brief Tensor storage re-optimization algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_OPTIMIZE
#define HAMMER_MATH_MULTIDIM_OPS_OPTIMIZE

#include <complex>

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class Optimize final {
        public:

            Optimize() noexcept = default;

            IContainer* operator()(VectorContainer& a);
            IContainer* operator()(SparseContainer& a);
            IContainer* operator()(OuterContainer& a);

            static IContainer* error(IContainer& /*a*/);
        };

        inline bool shouldBeSparse(size_t fill, size_t total) {
            return (fill * (sizeof(std::complex<double>) + sizeof(PositionType)) <
                    total * sizeof(std::complex<double>));
        }


    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
