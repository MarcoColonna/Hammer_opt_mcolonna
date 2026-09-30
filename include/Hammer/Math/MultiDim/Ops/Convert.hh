///
/// @file  Convert.hh
/// @brief Tensor storage type conversion algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_MULTIDIM_OPS_Convert
#define HAMMER_MATH_MULTIDIM_OPS_Convert

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer::MultiDimensional {

    class IContainer;
    class VectorContainer;
    class SparseContainer;
    class OuterContainer;

    namespace Ops {

        class Convert final {
        public:

            explicit Convert(bool destinationSparse = true) noexcept;

            static IContainer* error(IContainer& /*unused*/);

            IContainer* operator()(VectorContainer& a) const;
            IContainer* operator()(SparseContainer& a) const;
            IContainer* operator()(OuterContainer& a) const;
            IContainer* operator()(IContainer& a) const;

        private:

            static IContainer* toSparse(VectorContainer& a);
            static IContainer* toSparse(OuterContainer& a);

            static IContainer* toVector(SparseContainer& a);
            static IContainer* toVector(OuterContainer& a);

            bool _destIsSparse;
        };


    } // namespace Ops

} // namespace Hammer::MultiDimensional


#endif
