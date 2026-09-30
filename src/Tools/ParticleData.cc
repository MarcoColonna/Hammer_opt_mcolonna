///
/// @file  ParticleData.cc
/// @brief Hammer decay PDG code management
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <algorithm>
#include <iterator>

#include <boost/functional/hash.hpp>

#include "Hammer/Tools/ParticleData.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    ParticleData::~ParticleData() {
        _signatures.clear();
        _signatureIndex = 0ul;
    }

    HashId ParticleData::id() const {
        return _signatures[_signatureIndex].id;
    }

    HashId ParticleData::hadronicId() const {
        return _signatures[_signatureIndex].hadronicId;
    }

    const vector<double>& ParticleData::masses() const {
        return _signatures[_signatureIndex].masses;
    }

    bool ParticleData::checkValidSignature(PdgId parent, const vector<PdgId>& daughters,
                                           const vector<PdgId>& grandDaughters) {
        vector<PdgId> flippedDaughters;
        vector<PdgId> flippedGrandDaughters;
        flippedDaughters.reserve(daughters.size());
        flippedGrandDaughters.reserve(grandDaughters.size());
        transform(daughters.begin(), daughters.end(), back_inserter(flippedDaughters),
                  [parent](PdgId id) -> PdgId { return (parent > 0) ? id : flipSign(id); });
        transform(grandDaughters.begin(), grandDaughters.end(), back_inserter(flippedGrandDaughters),
                  [parent](PdgId id) -> PdgId { return (parent > 0) ? id : flipSign(id); });
        auto result1 = daughters;
        result1.reserve(daughters.size() + grandDaughters.size());
        result1.insert(result1.end(), grandDaughters.begin(), grandDaughters.end());
        auto result2 = combineDaughters(flippedDaughters, flippedGrandDaughters);
        transform(result2.begin(), result2.end(), result2.begin(),
                  [parent](PdgId id) -> PdgId { return (parent > 0) ? id : flipSign(id); });
        return equal(result1.begin(), result1.end(), result2.begin());
    }

    void ParticleData::addProcessSignature(PdgId parent, const vector<PdgId>& daughters,
                                           const vector<PdgId>& grandDaughters) {
        ASSERT_MSG(checkValidSignature(parent, daughters, grandDaughters),
                   "Particle order does not respect Hammer Convention in addProcessSignature. Shame on you!");
        Signature newsig{};
        newsig.parent = parent;
        newsig.daughters = combineDaughters(daughters, grandDaughters); // ordered by abs, for unique HashId
        newsig.calcIds(daughters.size());
        newsig.masses.clear();
        newsig.addMass(parent); // masses added in the order within addProcessSignature
        newsig.addMasses(daughters);
        newsig.addMasses(grandDaughters);
        _signatures.push_back(newsig);
    }

    void ParticleData::Signature::calcIds(size_t numDaughters) {
        id = processID(parent, daughters);
        vector<PdgId> hadDaughters;
        copy_n(daughters.begin(), numDaughters, back_inserter(hadDaughters));
        hadDaughters.erase(
            remove_if(hadDaughters.begin(), hadDaughters.end(), [](PdgId idx) { return (abs(idx) < 100); }),
            hadDaughters.end());
        hadronicId = processID(parent, hadDaughters);
    }

    void ParticleData::Signature::addMass(PdgId idx) {
        const PID& pdg = PID::instance();
        masses.push_back(pdg.getMass(idx));
    }

    void ParticleData::Signature::addMasses(const vector<PdgId>& ids) {
        const PID& pdg = PID::instance();
        for (auto elem : ids) {
            masses.push_back(pdg.getMass(elem));
        }
    }

    bool ParticleData::setSignatureIndex(size_t idx) {
        if (idx > _signatures.size()) {
            _signatureIndex = 0ul;
            return false;
        }
        _signatureIndex = idx;
        return true;
    }

    size_t ParticleData::numSignatures() const {
        return _signatures.size();
    }

} // namespace Hammer
