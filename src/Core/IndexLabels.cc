///
/// @file  IndexLabels.cc
/// @brief Tensor indices label function implementations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include <cstdlib>

#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;

namespace Hammer {

    bool isSpinIndex(IndexLabel val) {
        return (abs(val) <= SPIN_INDEX_END);
    }

    bool isFFIndex(IndexLabel val) {
        auto posval = abs(val);
        return ((posval >= FF_INDEX_START) && (posval <= FF_VAR_INDEX_END));
    }

    bool isWCIndex(IndexLabel val) {
        auto posval = abs(val);
        return ((posval >= WC_INDEX_START) && (posval <= WC_INDEX_END));
    }

    bool isSpecializedIndex(IndexLabel val) {
        return (abs(val) > SPECIALIZATION_PAD);
    }

    pair<uint32_t, int> getSpecializedLabelComponents(IndexLabel val) {
        auto tmp = div(val, SPECIALIZATION_PAD);
        return make_pair(static_cast<uint32_t>(abs(tmp.quot)), static_cast<int>(tmp.rem));
    }

    IndexLabel specializeLabel(IndexLabel val, uint32_t pad) {
        return static_cast<IndexLabel>(val + ((val > NONE) ? 1 : -1) * SPECIALIZATION_PAD * pad);
    }

    uint32_t generateSpecializedPad(const std::string& id) {
        return Hash::murmur3_32(reinterpret_cast<const uint8_t*>(id.c_str()), id.size(), 31ul);
    }


} // namespace Hammer
