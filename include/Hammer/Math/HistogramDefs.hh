#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
///
/// @file  HistogramDefs.hh
/// @brief Hammer histogram class template methods definitions
///
#pragma clang diagnostic pop

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include "Hammer/Tools/Logging.hh"

namespace Hammer {

    template <typename... Args>
    void Histogram::fillBin(const Tensor& value, Args... rest) {
        try {
            element({static_cast<IndexType>(rest)...}).addWeight(value);
        } catch (RangeError& e) {
            MSG_ERROR(std::string("Error setting element: ") + e.what() + std::string(" Element not set."));
        }
    }

    template <typename... Args>
    const Tensor& Histogram::getWeight(Args... rest) const {
        return element({static_cast<IndexType>(rest)...}).getWeight();
    }

    template <typename... Args>
    const Tensor& Histogram::getSquaredWeight(Args... rest) const {
        return element({static_cast<IndexType>(rest)...}).getSquaredWeight();
    }

    template <typename... Args>
    size_t Histogram::getNumEvents(Args... rest) const {
        return element({static_cast<IndexType>(rest)...}).getNumEvents();
    }

} // namespace Hammer
