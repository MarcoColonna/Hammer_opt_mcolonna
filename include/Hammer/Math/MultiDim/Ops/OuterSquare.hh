///
/// @file  OuterSquare.hh
/// @brief Tensor outer square algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_OUTERSQUARE
#define HAMMER_MATH_MULTIDIM_OPS_OUTERSQUARE

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class OuterSquare final {
        public:

            OuterSquare() noexcept = default;

            IContainer* operator()(OuterContainer& a);

            IContainer* operator()(IContainer& a);

            static IContainer* error(IContainer& /*unused*/);
        };

    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
