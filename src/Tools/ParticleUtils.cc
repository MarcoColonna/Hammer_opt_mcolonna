///
/// @file  ParticleUtils.cc
/// @brief PDG codes to UID functions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <algorithm>

#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    std::vector<PdgId> flipSigns(const std::vector<PdgId>& list) {
        vector<PdgId> res = list;
        for (auto& elem : res) {
            // check to ensure particle is not self conjugate particle: isMeson return false for -ve pid self conj
            // meson.
            if (Hammer::PID::isMeson(elem) == Hammer::PID::isMeson(-elem) && elem != PID::PHOTON) {
                elem *= -1;
            }
        }
        return res;
    }

    PdgId flipSign(const PdgId& id) {
        PdgId res = id;
        // check to ensure particle is not self conjugate meson: isMeson return false for -ve pid self conj meson.
        if (Hammer::PID::isMeson(res) == Hammer::PID::isMeson(-res) && res != PID::PHOTON) {
            res *= -1;
        }
        return res;
    }

    std::vector<PdgId> combineDaughters(const std::vector<PdgId>& daughters, const std::vector<PdgId>& subDaughters) {
        vector<PdgId> res = daughters;
        sort(res.begin(), res.end(), &pdgSorter);
        if (!subDaughters.empty()) {
            vector<PdgId> tmp = subDaughters;
            sort(tmp.begin(), tmp.end(), &pdgSorter);
            res.insert(res.end(), tmp.begin(), tmp.end());
        }
        return res;
    }

    HashId processID(PdgId parent, const std::vector<PdgId>& allDaughters) {
        HashId seed = 0ul;
        combine_hash(seed, static_cast<int>(parent));
        for (auto elem : allDaughters) {
            combine_hash(seed, elem);
        }
        return seed;
    }

    HashId combineProcessIDs(const std::set<HashId>& allIds) {
        HashId seed = 0ul;
        for (auto elem : allIds) {
            combine_hash(seed, elem);
        }
        return seed;
    }


    bool pdgSorter(PdgId a, PdgId b) {
        if (abs(static_cast<int>(a)) != abs(static_cast<int>(b))) {
            return (abs(static_cast<int>(a)) > abs(static_cast<int>(b)));
        }
        return (a > b);
    }

    bool particlesByPdg(const std::function<PdgId(const Particle&)>& pdgGetter, const Particle& a, const Particle& b) {
        PdgId aId = pdgGetter(a);
        PdgId bId = pdgGetter(b);
        return pdgSorter(aId, bId);
    }

} // namespace Hammer
