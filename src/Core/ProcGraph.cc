///
/// @file  ProcGraph.cc
/// @brief Container class for process tree structure and its amplitudes associations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/ProcGraph.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Process.hh"

using namespace std;

namespace Hammer {

    AmplEntry AmplTriplet::resolvePurePS(bool parentPS, bool daughterPS) const {
        int switchVal = (parent == nullptr ? 0 : 1) * 8 + (daughter == nullptr ? 0 : 1) * 4 + (parentPS ? 1 : 0) * 2 +
                        (daughterPS ? 1 : 0);
        if (switchVal % 4 == 0) { // no PS
            return AmplEntry{edge, AmplType::FULLEDGE};
        }
        if (switchVal % 4 == 3) { // pure PS edge
            return AmplEntry{nullptr, AmplType::FULLEDGE};
        }
        if ((abs(switchVal - 3) == 2)) { // daughter declared PS, but not enforced
            MSG_ERROR("Daughter of an edge declared PS with undeclared parent. You can't touch this!");
            return AmplEntry{edge, AmplType::FULLEDGE};
        }
        if ((abs(switchVal - 6) == 4)) { // parent declared PS, but not enforced
            MSG_ERROR("Parent of an edge declared PS with undeclared daughter. You can't touch this!");
            return AmplEntry{edge, AmplType::FULLEDGE};
        }
        if ((abs(switchVal - 11) == 2)) { // daughter declared PS; only parent partial edge returned
            return AmplEntry{parent, AmplType::PARENTEDGE};
        }
        if ((abs(switchVal - 10) == 4)) { // parent declared PS; only daughter partial edge returned
            return AmplEntry{daughter, AmplType::DAUGHTEREDGE};
        }
        MSG_ERROR("I can't count!");
        return AmplEntry{edge, AmplType::FULLEDGE};
    }

    Log& AmplTriplet::getLog() {
        return Log::getLog("Hammer.AmplTriplet");
    }

    void ProcGraph::initialize(const DictionaryManager* dictionaries, const Process* inputs) {
        _dictionaries = dictionaries;
        _inputs = inputs;
        makeDepthMap(_inputs->getFirstVertex(), 0ul);
        makeVertexEdgeTrees();
        selectEdgeVertices();
    }

    const ProcGraph::SelAmplitudeDict& ProcGraph::selectedAmplitudes() const {
        return _selectedAmplitudes;
    }
    const UniqueParticleIndices& ProcGraph::assignedVertices() const {
        return _assignedVertices;
    }

    void ProcGraph::makeDepthMap(ParticleIndex parent, size_t seed) {
        auto it = _inputs->find(parent);
        if (it == _inputs->end()) {
            return;
        }
        _depthMap.insert({parent, seed});
        for (const auto& daughter : it->second) {
            makeDepthMap(daughter, seed + 1ul);
        }
    }

    void ProcGraph::makeVertexEdgeTrees() {
        _implementedVertices.clear();
        _implementedEdges.clear();
        // std::map<ParticleIndex, size_t> depthMap;
        // depthMap.clear();
        ////starts from topmost vertex in the processTree
        // depthMap.insert({_inputs->getFirstVertex(), 0ul});
        for (const auto& elem : *_inputs) {
            vector<PdgId> daughterPdgs;
            size_t parWeight = 0ul;
            daughterPdgs.reserve(elem.second.size());
            for (auto elem2 : elem.second) {
                daughterPdgs.push_back(_inputs->getParticle(elem2).pdgId());
            }
            // check and fill implemented vertex, whether to be evaluated or purePS
            PdgId pdgParent = _inputs->getParticle(elem.first).pdgId();
            auto* ampl = _dictionaries->providers().getAmplitude(pdgParent, daughterPdgs);
            auto purePS = _dictionaries->purePSDefs().isPurePhaseSpace(pdgParent, daughterPdgs);
            if (ampl != nullptr) {
                VertexEntry amplTupl{elem.first, ampl, purePS};
                _implementedVertices.push_back(amplTupl);
                parWeight = 1ul;
            }
            // check and fill implemented edge, whether to be evaluated or purePS
            for (auto elem2 : elem.second) {
                auto it = _inputs->find(elem2);
                if (it != _inputs->end()) {
                    vector<PdgId> granddaughterPdgs;
                    // depthMap.insert({it->first, depthMap.find(elem.first)->second + 1ul});
                    for (auto elem3 : it->second) {
                        PdgId pdgCodeGrandDaughter = _inputs->getParticle(elem3).pdgId();
                        granddaughterPdgs.push_back(pdgCodeGrandDaughter);
                    }
                    auto* edgeAmpl =
                        _dictionaries->providers().getAmplitude(pdgParent, daughterPdgs, granddaughterPdgs);
                    if (edgeAmpl != nullptr) {
                        size_t edgeWeight = parWeight + 1ul;
                        PdgId pdgCodeDaughter = _inputs->getParticle(elem2).pdgId();
                        auto* daughterAmpl =
                            _dictionaries->providers().getAmplitude(pdgCodeDaughter, granddaughterPdgs);
                        auto daughterPurePS =
                            _dictionaries->purePSDefs().isPurePhaseSpace(pdgCodeDaughter, granddaughterPdgs);
                        if (daughterAmpl != nullptr) {
                            edgeWeight += 1ul;
                        }
                        AmplTriplet amplTriplet{ampl, daughterAmpl, edgeAmpl};
                        EdgeEntry edgeTupl{elem.first, elem2, amplTriplet, purePS, daughterPurePS};
                        _implementedEdges[edgeWeight][_depthMap.find(elem2)->second].push_back(edgeTupl);
                    }
                }
            }
        }
    }

    bool ProcGraph::isAssignedVertex(ParticleIndex index) const {
        return _assignedVertices.find(index) != _assignedVertices.end();
    }

    void ProcGraph::selectEdgeVertices() {
        _assignedVertices.clear();
        _selectedAmplitudes.clear();
        for (auto& edgeWeight : _implementedEdges) {
            for (auto& depthLevel : edgeWeight.second) { // this map is ordered largest to smallest
                for (auto& amplTuple : depthLevel.second) {
                    if (!isAssignedVertex(amplTuple.parentIdx) && !isAssignedVertex(amplTuple.daughterIdx)) {
                        auto amplTriplet = amplTuple.amplitudes; // this is a struct of 3 AmplitudeBase ptrs
                        AmplEntry numAmpl = amplTriplet.resolvePurePS(amplTuple.parentPSFlags.numerator,
                                                                      amplTuple.daughterPSFlags.numerator);
                        AmplEntry denAmpl = amplTriplet.resolvePurePS(amplTuple.parentPSFlags.denominator,
                                                                      amplTuple.daughterPSFlags.denominator);
                        SelectedAmplEntry tmpampl{NumDenPair<AmplEntry>{numAmpl, denAmpl}, amplTuple.daughterIdx};
                        _selectedAmplitudes.insert({amplTuple.parentIdx, tmpampl});
                        _assignedVertices.insert(amplTuple.parentIdx);
                        _assignedVertices.insert(amplTuple.daughterIdx);
                    }
                }
            }
        }
        for (auto& amplTuple : _implementedVertices) {
            if (!isAssignedVertex(amplTuple.parentIdx)) {
                AmplEntry numAmpl{amplTuple.parentPSFlags.numerator ? nullptr : amplTuple.amplitude, AmplType::VERTEX};
                AmplEntry denAmpl{amplTuple.parentPSFlags.denominator ? nullptr : amplTuple.amplitude,
                                  AmplType::VERTEX};
                SelectedAmplEntry tmpampl{NumDenPair<AmplEntry>{numAmpl, denAmpl}, _inputs->getFirstVertex()};
                _selectedAmplitudes.insert({amplTuple.parentIdx, tmpampl});
                _assignedVertices.insert(amplTuple.parentIdx);
            }
        }
    }

    tuple<Particle, ParticleList, ParticleList> ProcGraph::getParticleVectors(AmplType ampltype, ParticleIndex parent,
                                                                              ParticleIndex daughter) const {
        if (!_inputs->isParent(parent) || !_inputs->isParent(daughter)) {
            MSG_ERROR("Something really stinky just happened in CalcAmplitudes. Can you smell it?");
        }
        switch (ampltype) {
        case AmplType::VERTEX: {
            return tuple<Particle, ParticleList, ParticleList>{
                _inputs->getParticle(parent), _inputs->getDaughters(parent, true), _inputs->getSiblings(parent, true)};
        }
        case AmplType::FULLEDGE: {
            auto daughters = _inputs->getDaughters(parent, true);
            auto tmp = _inputs->getDaughters(daughter, true);
            daughters.insert(daughters.end(), tmp.begin(), tmp.end());
            return tuple<Particle, ParticleList, ParticleList>{_inputs->getParticle(parent), daughters,
                                                               _inputs->getSiblings(parent, true)};
        }
        case AmplType::PARENTEDGE: {
            return tuple<Particle, ParticleList, ParticleList>{
                _inputs->getParticle(parent), _inputs->getDaughters(parent, true), _inputs->getSiblings(parent, true)};
        }
        case AmplType::DAUGHTEREDGE: {
            return tuple<Particle, ParticleList, ParticleList>{_inputs->getParticle(daughter),
                                                               _inputs->getDaughters(daughter, true),
                                                               _inputs->getSiblings(daughter, true)};
        }
        }
        MSG_ERROR("Something really stinky just happened in CalcAmplitudes. Can you smell it?");
        return tuple<Particle, ParticleList, ParticleList>{
            _inputs->getParticle(parent), _inputs->getDaughters(parent, true), _inputs->getSiblings(parent, true)};
    }

    double ProcGraph::getMassFactor(AmplEntry ampl, ParticleIndex parent, ParticleIndex daughter) const {
        if (!_inputs->isParent(parent) || !_inputs->isParent(daughter)) {
            MSG_ERROR("Something really stinky just happened in getMassFactor. Can you smell it?");
        }
        auto pmass = _inputs->getParticle(parent).p().mass();
        auto dsize = static_cast<double>(_inputs->getDaughters(parent).size());
        switch (ampl.type) {
        case AmplType::VERTEX: {
            if (ampl.ptr == nullptr) { // check if PS vertex
                return pow(pmass, 6. - 2. * dsize);
            }
            break;
        }
        case AmplType::FULLEDGE: {
            if (ampl.ptr == nullptr) { // check if PS edge, i.e. both vertices are PS
                auto dmass = _inputs->getParticle(daughter).p().mass();
                auto gdsize = static_cast<double>(_inputs->getDaughters(daughter).size());
                return pow(pmass, 6. - 2. * dsize) * pow(dmass, 6. - 2. * gdsize);
            }
            break;
        }
        case AmplType::PARENTEDGE: {
            auto dmass = _inputs->getParticle(daughter).p().mass();
            auto gdsize = static_cast<double>(_inputs->getDaughters(daughter).size());
            return pow(dmass, 6. - 2. * gdsize);
        }
        case AmplType::DAUGHTEREDGE: {
            return pow(pmass, 6. - 2. * dsize);
        }
        }
        return 1.;
    }

    double ProcGraph::getMultFactor(AmplEntry ampl, AmplEntry ancestorAmpl, bool isAncestorGrandparent) {
        auto val = static_cast<double>(ampl.ptr->multiplicityFactor());
        if (ampl.type == AmplType::DAUGHTEREDGE) { // No matter the antecedent in the process tree, inside the edge, the
                                                   // daughter has PS parent.
            return 1. / val;
        }
        if (isAncestorGrandparent) { // Amplitude is connected to the daughter vertex of an edge
            if (ancestorAmpl.ptr == nullptr || ancestorAmpl.type == AmplType::PARENTEDGE) {
                return 1. / val;
            }
        } else { // Amplitude is connected to the parent vertex of a daughter edge or just a simple vertex
            if (ancestorAmpl.ptr == nullptr || ancestorAmpl.type == AmplType::DAUGHTEREDGE) {
                return 1. / val;
            }
        }
        return 1.;
    }

    Log& ProcGraph::getLog() {
        return Log::getLog("Hammer.ProcGraph");
    }

} // namespace Hammer
