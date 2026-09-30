///
/// @file  Process.cc
/// @brief Hammer process class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/functional/hash.hpp>

#include "Hammer/Process.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;

namespace Hammer {

    Process::Process() = default;

    // Process::~Process() noexcept {
    //     _processTree.clear();
    //     _particles.clear();
    //     _parentDict.clear();
    //     _idTree.clear();
    //     _fullId.clear();
    //     _erasedGammas.clear();
    // }

    Process::const_iterator Process::begin() const {
        return _processTree.begin();
    }

    Process::const_iterator Process::end() const {
        return _processTree.end();
    }

    Process::const_iterator Process::find(ParticleIndex particle) const {
        return _processTree.find(particle);
    }

    ParticleIndex Process::addParticle(const Particle& p) {
        _initialized = false;
        ParticleIndex tmp = _particles.size();
        _particles.push_back(p);
        return tmp;
    }

    void Process::addVertex(ParticleIndex parent, const ParticleIndices& daughters) {
        _initialized = false;
        bool can_add = parent < _particles.size();
        for (auto elem : daughters) {
            can_add &= elem < _particles.size();
        }
        if (can_add) {
            _processTree[parent] = daughters;
        } else {
            MSG_ERROR("One or more particle index is invalid: vertex NOT added.");
        }
    }

    void Process::removeVertex(ParticleIndex parent, bool prune) {
        _initialized = false;
        auto it = _processTree.find(parent);
        if (it != _processTree.end()) {
            auto daughters = daughtersId(it);
            if (prune) {
                for (auto elem : daughters) {
                    removeVertex(elem, true);
                }
            } else {
                if (parent == 0) {
                    MSG_ERROR("You can't collapse after removing the tree head!");
                    return;
                }
                auto it2 = _processTree.find(getParentId(parent));
                if (it2 != _processTree.end()) {
                    it2->second.insert(it2->second.end(), daughters.begin(), daughters.end());
                }
            }
            _processTree.erase(parent);
        }
    }

    const ParticleIndices& Process::getDaughtersIds(ParticleIndex parent) const {
        return getOrThrow(_processTree, parent, RangeError("Index not found"));
    }

    ParticleIndices Process::getSiblingsIds(ParticleIndex particle) const {
        ParticleIndices temps;
        ParticleIndex parent = getParentId(particle);
        if (parent >= _particles.size()) {
            return temps;
        }
        auto res = getDaughtersIds(parent);
        copy_if(res.begin(), res.end(), back_inserter(temps),
                [particle](const ParticleIndex& elem) -> bool { return elem != particle; });
        return temps;
    }

    ParticleIndex Process::getParentId(ParticleIndex daughter) const {
        if (_initialized) {
            return getOrDefault(_parentDict, daughter, numeric_limits<size_t>::max());
        }
        for (const auto& elem : _processTree) {
            if (::find(elem.second.begin(), elem.second.end(), daughter) != elem.second.end()) {
                return elem.first;
            }
        }
        return numeric_limits<size_t>::max();
    }

    const Particle& Process::getParticle(ParticleIndex id) const {
        if (id < _particles.size()) {
            return _particles[id];
        }
        throw Error("Index out of range. You're gonna need a bigger boat, or a smaller shark.");
    }

    ParticleList Process::sortByParent(ParticleIndex parent, ParticleList parts) const {
        PdgId parentIdVal = _particles[parent].pdgId();
        function<PdgId(const Particle&)> getter = [parentIdVal](const Particle& a) -> PdgId {
            return (parentIdVal > 0) ? a.pdgId() : -a.pdgId();
        };
        auto sorter = [&getter](auto&& PH1, auto&& PH2) {
            return particlesByPdg(getter, std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2));
        };
        sort(parts.begin(), parts.end(), sorter);
        return parts;
    }

    ParticleList Process::getDaughters(ParticleIndex parent, bool sorted) const {
        ParticleList temps;
        auto res = getDaughtersIds(parent);
        temps.reserve(res.size());
        transform(res.begin(), res.end(), back_inserter(temps),
                  [this](const ParticleIndex& elem) -> const Particle& { return _particles[elem]; });
        if (sorted) {
            return sortByParent(parent, temps);
        }
        return temps;
    }

    const Particle& Process::getParent(ParticleIndex daughter) const {
        ParticleIndex tmp = getParentId(daughter);
        return _particles[tmp];
    }

    ParticleList Process::getSiblings(ParticleIndex particle, bool sorted) const {
        ParticleIndex parent = getParentId(particle);
        if (parent >= _particles.size()) {
            return {};
        }
        ParticleList temps;
        auto res = getDaughtersIds(parent);
        temps.reserve(res.size() - 1);
        for (auto& elem : res) {
            if (elem != particle) {
                temps.push_back(_particles[elem]);
            }
        }
        if (sorted) {
            return sortByParent(parent, temps);
        }
        return temps;
    }

    bool Process::isParent(ParticleIndex particle) const {
        return _processTree.find(particle) != _processTree.end();
    }

    pair<Particle, ParticleList> Process::getParticlesByVertex(PdgId parent, const vector<PdgId>& daughters) const {
        HashId vertexId = processID(parent, combineDaughters(daughters));
        auto it = _idTree.find(vertexId);
        if (it != _idTree.end()) {
            return {_particles[it->second], getDaughters(it->second, true)};
        }
        vertexId = processID(flipSign(parent), flipSigns(combineDaughters(daughters)));
        it = _idTree.find(vertexId);
        if (it != _idTree.end()) {
            return {_particles[it->second], getDaughters(it->second, true)};
        }
        return {};
    }

    pair<Particle, ParticleList> Process::getParticlesByVertex(const string& vertex) const {
        const PID& pdg = PID::instance();
        auto temp = pdg.expandToValidVertexUIDs(vertex);
        for (auto& elem : temp) {
            auto it = _idTree.find(elem);
            if (it != _idTree.end()) {
                return {_particles[it->second], getDaughters(it->second, true)};
            }
        }
        return {};
    }

    HashId Process::getVertexId(const string& vertex) const {
        const PID& pdg = PID::instance();
        auto temp = pdg.expandToValidVertexUIDs(vertex);
        auto it =
            find_if(temp.begin(), temp.end(), [&](HashId val) -> bool { return _fullId.find(val) != _fullId.end(); });
        if (it == temp.end()) {
            return 0;
        }
        return *it;
        // for (auto& elem : temp) {
        //     auto it = _fullId.find(elem);
        //     if (it != _fullId.end()){
        //         return *it;
        //     }
        // }
        // return 0;
    }

    size_t Process::numParticles(bool withoutPhotons) const {
        return _particles.size() - (withoutPhotons ? _erasedGammas.size() : 0ul);
    }

    ParticleIndex Process::getFirstVertex() const {
        if (!_initialized) {
            return findFirstVertex();
        }
        return _firstVertexIdx;
    }

    ParticleIndex Process::findFirstVertex() const {
        UniqueParticleIndices used = _erasedGammas;
        for (const auto& elem : _processTree) {
            for (const auto& elem2 : elem.second) {
                used.insert(elem2);
            }
        }
        for (ParticleIndex i = 0; i < _particles.size(); ++i) {
            if (used.find(i) == used.end()) {
                return i;
            }
        }
        return numeric_limits<size_t>::max();
    }

    void Process::cacheParticleDependencies() {
        UniqueParticleIndices used = _erasedGammas;
        _parentDict.reserve(_particles.size() - 1);
        for (auto& elem : _processTree) {
            for (auto& elem2 : elem.second) {
                _parentDict.insert({elem2, elem.first});
                used.insert(elem2);
            }
        }
        for (ParticleIndex i = 0; i < _particles.size(); ++i) {
            if (used.find(i) == used.end()) {
                _firstVertexIdx = i;
                break;
            }
        }
    }


    bool Process::initialize() {
        if (_particles.empty()) {
            return false;
        }
        if (!_initialized) {
            pruneSoftPhotons();
            calcSignatures();
            cacheParticleDependencies();
            _initialized = true;
        }
        return true;
        /// @todo what else needs to be done here?
        // The 'else'. We'll never get that two hours back.
    }

    void Process::pruneSoftPhotons() {
        auto firstvertex = findFirstVertex();
        for (auto& elem : _processTree) {
            auto& daughters = elem.second;
            UniqueParticleIndices erasedgammas = _erasedGammas; // includes gammas deleted at other vertices
            while (daughters.size() > 2 && find_if(daughters.begin(), daughters.end(), [&](ParticleIndex& idx) -> bool {
                                               return (_particles[idx].pdgId() == PID::PHOTON);
                                           }) != daughters.end()) {
                // get softest photon
                map<double, pair<ParticleIndex, Particle>> gammas;
                for (auto idx : daughters) {
                    if (_particles[idx].pdgId() == PID::PHOTON) {
                        gammas.insert({_particles[idx].p().E(), make_pair(idx, _particles[idx])});
                    }
                }
                pair<ParticleIndex, Particle> gamma = gammas.begin()->second;
                // get nearest charged particle
                map<double, pair<ParticleIndex, Particle>, std::greater<>> chps;
                for (auto idx : daughters) {
                    if (PID::getThreeCharge(_particles[idx].pdgId()) != 0) {
                        const auto& mom1 = _particles[idx].p();
                        const auto& mom2 = gamma.second.p();
                        chps.insert({costheta(mom1.pVec(), mom2.pVec()), make_pair(idx, _particles[idx])});
                    }
                }
                pair<ParticleIndex, Particle> nchp = chps.begin()->second;
                // get 'Y' constituents
                set<ParticleIndex> chTree{nchp.first};
                ParticleIndices siblings = getSiblingsIds(nchp.first);
                siblings.erase(remove(siblings.begin(), siblings.end(), gamma.first), siblings.end());
                set<ParticleIndex> YTree(siblings.begin(), siblings.end());
                // Momenta
                FourMomentum pY;
                for (auto elemY : YTree) {
                    pY += _particles[elemY].p();
                }
                FourMomentum pCh = nchp.second.p();
                // Get sub tree members
                for (ParticleIndex idx = 0; idx < _particles.size(); ++idx) {
                    auto ix = idx;
                    while ((ix != firstvertex) && (::find(daughters.begin(), daughters.end(), ix) == daughters.end()) &&
                           (erasedgammas.find(ix) == erasedgammas.end())) {
                        ix = getParentId(ix);
                    }
                    if (ix == nchp.first) {
                        chTree.insert(idx);
                    }
                    if (::find(siblings.begin(), siblings.end(), ix) != siblings.end() && ix != gamma.first) {
                        YTree.insert(idx);
                    }
                }
                // Now apply boosts. P rest frame objects
                FourMomentum pP = _particles[elem.first].p();
                double pYCh = 0.5 * pP.mass() *
                              sqrt((1. - pow((pCh.mass() + pY.mass()), 2.) / pP.mass2()) *
                                   (1. - pow((pCh.mass() - pY.mass()), 2.) / pP.mass2()));
                double EY = (pP.mass2() - pCh.mass2() + pY.mass2()) / (2. * pP.mass());
                double ECh = (pP.mass2() + pCh.mass2() - pY.mass2()) / (2. * pP.mass());
                // Boost to P frame
                auto boostP = pY + pCh;
                pY.boostToRestFrameOf(boostP);
                pCh.boostToRestFrameOf(boostP);
                // Construct separate boosts on Ch and Y to restore momentum conservation
                double betagammaCh = (pYCh * pCh.E() - ECh * pCh.p()) / pCh.mass2();
                double betagammaY = (pYCh * pY.E() - EY * pY.p()) / pY.mass2();
                double betaCh = (betagammaCh > 0.0) ? 1. / sqrt(1. + 1. / pow(betagammaCh, 2.))
                                                    : -1. / sqrt(1. + 1. / pow(betagammaCh, 2.));
                double betaY = (betagammaY > 0.0) ? 1. / sqrt(1. + 1. / pow(betagammaY, 2.))
                                                  : -1. / sqrt(1. + 1. / pow(betagammaY, 2.));
                FourMomentum boostCh{1., betaCh * pCh.px() / pCh.p(), betaCh * pCh.py() / pCh.p(),
                                     betaCh * pCh.pz() / pCh.p()};
                FourMomentum boostY{1., betaY * pY.px() / pY.p(), betaY * pY.py() / pY.p(), betaY * pY.pz() / pY.p()};
                for (auto elemCh : chTree) {
                    _particles[elemCh]
                        .p()
                        .boostToRestFrameOf(boostP)
                        .boostFromRestFrameOf(boostCh)
                        .boostFromRestFrameOf(pP);
                }
                for (auto elemY : YTree) {
                    _particles[elemY].p().boostToRestFrameOf(boostP).boostFromRestFrameOf(boostY).boostFromRestFrameOf(
                        pP);
                }
                daughters.erase(remove(daughters.begin(), daughters.end(), gamma.first), daughters.end());
                erasedgammas.insert(gamma.first);
            }
            _erasedGammas.insert(erasedgammas.begin(), erasedgammas.end());
        }
    }


    void Process::calcSignatures() {
        _fullId.clear();
        _idTree.clear();
        for (auto& elem : _processTree) {
            vector<PdgId> daughterPdgs;
            transform(elem.second.begin(), elem.second.end(), back_inserter(daughterPdgs),
                      [this](const ParticleIndex& elem2) -> PdgId { return _particles[elem2].pdgId(); });
            if (::find(daughterPdgs.begin(), daughterPdgs.end(), _particles[elem.first].pdgId()) !=
                daughterPdgs.end()) {
                continue;
            }
            HashId decayId = processID(_particles[elem.first].pdgId(), combineDaughters(daughterPdgs));
            _fullId.insert(decayId);
            _idTree.insert({decayId, elem.first});
        }
        _hashId = combineProcessIDs(_fullId);
    }

    HashId Process::getId() const {
        return _hashId;
    }

    const set<HashId>& Process::fullId() const {
        return _fullId;
    }

    void Process::disable() {
        _hashId = 0;
    }


    void Process::write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBProcIDs>* msg) const {
        vector<uint64_t> idfull;
        for (auto elem : _fullId) {
            idfull.push_back(elem);
        }
        auto serialfull = msgwriter->CreateVector(idfull);
        Serial::FBProcIDsBuilder serialprocess{*msgwriter};
        serialprocess.add_id(_hashId);
        serialprocess.add_idfull(serialfull);
        *msg = serialprocess.Finish();
    }

    void Process::read(const Serial::FBProcIDs* msgreader) {
        if (msgreader != nullptr) {
            _hashId = msgreader->id();
            const auto* idfull = msgreader->idfull();
            if (idfull == nullptr) {
                MSG_ERROR("Malformed read: you need a song and sing it in the morning!");
                return;
            }
            _fullId.clear();
            for (unsigned int i = 0; i < idfull->size(); ++i) {
                _fullId.insert(idfull->Get(i));
            }
        }
    }

    void Process::defineSettings() {
        setPath("Process");
    }

    Log& Process::getLog() {
        return Log::getLog("Hammer.Process");
    }

} // namespace Hammer
