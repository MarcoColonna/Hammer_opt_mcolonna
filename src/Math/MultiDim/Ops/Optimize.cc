///
/// @file  Optimize.cc
/// @brief Tensor storage re-optimization algorithm
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Math/MultiDim/Ops/Optimize.hh"
#include "Hammer/Math/MultiDim/Ops/Convert.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Math/MultiDim/OuterContainer.hh"

using namespace std;

namespace Hammer::MultiDimensional {

    using VTensor = VectorContainer;
    using STensor = SparseContainer;
    using OTensor = OuterContainer;
    using Base = IContainer;

    namespace Ops {

        Base* Optimize::operator()(STensor& a) {
            if (shouldBeSparse(a.dataSize() / a.entrySize(), a.numValues())) {
                return (Convert{false})(a);
            }
            return static_cast<Base*>(&a);
        }

        Base* Optimize::operator()(VTensor& a) {
            if (shouldBeSparse(a.dataSize() / a.entrySize(), a.numValues())) {
                return (Convert{true})(a);
            }
            return static_cast<Base*>(&a);
        }

        Base* Optimize::operator()(OTensor& a) {
            // first check whether we need to collapse it to sparse or vector tensors
            if (a.shouldBeEvaluated()) {
                auto* sparse = dynamic_cast<STensor*>((Convert{true})(a));
                return (*this)(*sparse);
            }
            // otherwise loop over constituents and optimize those
            auto ptrs = a.getUniquePtrs();
            for (auto* elem : ptrs) {
                Base* res = nullptr;
                auto* sparse = dynamic_cast<STensor*>(elem);
                /// @todo improve here with auto dispatch if more tensors will be added
                if (sparse != nullptr) {
                    res = (*this)(*sparse);
                } else {
                    auto* vec = dynamic_cast<VTensor*>(elem);
                    if (vec != nullptr) {
                        res = (*this)(*vec);
                    }
                }
                if (res != nullptr && res != elem) {
                    a.swapElement(elem, TensorData{res});
                }
            }
            return static_cast<Base*>(&a);
        }

        Base* Optimize::error(Base& a) {
            return &a;
        }

    } // namespace Ops


} // namespace Hammer::MultiDimensional
