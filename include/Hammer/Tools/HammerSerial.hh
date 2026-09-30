///
/// @file  HammerSerial.hh
/// @brief Serialization related typedefs and includes
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_HammerSerial_HH
#define HAMMER_HammerSerial_HH

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wcovered-switch-default"
#pragma clang diagnostic ignored "-Wextra-semi"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wundef"
#pragma clang diagnostic ignored "-Wreserved-id-macro"

#include "Hammer/Tools/HammerSerial_generated.h"

#pragma clang diagnostic pop

namespace Hammer::Serial {

    class DetachedBuffers {
    public:

        DetachedBuffers() = default;
        DetachedBuffers(const DetachedBuffers&) = delete;
        DetachedBuffers& operator=(const DetachedBuffers&) = delete;
        DetachedBuffers(DetachedBuffers&&) = default;
        DetachedBuffers& operator=(DetachedBuffers&&) = default;
        ~DetachedBuffers() = default;

        void add(flatbuffers::DetachedBuffer&& elem, char kind);

        [[nodiscard]] size_t size() const;

        void clear();

        [[nodiscard]] const flatbuffers::DetachedBuffer& buffer(size_t pos) const;
        [[nodiscard]] char bufferType(size_t pos) const;

    private:

        std::vector<flatbuffers::DetachedBuffer> _data;
        std::vector<char> _types;
    };

} // namespace Hammer::Serial

#endif
