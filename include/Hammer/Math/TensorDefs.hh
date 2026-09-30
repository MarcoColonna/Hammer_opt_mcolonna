#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
///
/// @file  TensorDefs.hh
/// @brief Hammer tensor class template methods definitions
///
#pragma clang diagnostic pop

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include "Hammer/Tools/Logging.hh"
#include "Hammer/Exceptions.hh"

namespace Hammer {

    template <typename... Args>
    void Tensor::setValue(std::complex<double> value, Args... rest) {
        try {
            element({static_cast<IndexType>(rest)...}) = value;
        } catch (RangeError& e) {
            MSG_ERROR(std::string("Error setting element: ") + e.what() + std::string(" Element not set."));
        }
    }

    template <typename... Args>
    std::complex<double> Tensor::value(Args... rest) const {
        return element({static_cast<IndexType>(rest)...});
    }


} // namespace Hammer
