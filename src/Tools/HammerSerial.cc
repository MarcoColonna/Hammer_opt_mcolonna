///
/// @file  HammerSerial.cc
/// @brief Serialization related typedefs and includes
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Tools/HammerSerial.hh"

using namespace std;

namespace Hammer::Serial {

    void DetachedBuffers::add(flatbuffers::DetachedBuffer&& elem, char kind) {
        _data.push_back(std::move(elem));
        _types.push_back(kind);
    }

    size_t DetachedBuffers::size() const {
        return _data.size();
    }

    void DetachedBuffers::clear() {
        _data.clear();
        _types.clear();
    }

    const flatbuffers::DetachedBuffer& DetachedBuffers::buffer(size_t pos) const {
        return _data.at(pos);
    }

    char DetachedBuffers::bufferType(size_t pos) const {
        return _types.at(pos);
    }

} // namespace Hammer::Serial
