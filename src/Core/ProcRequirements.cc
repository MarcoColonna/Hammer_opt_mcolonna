///
/// @file  ProcRequirements.cc
/// @brief Container class for required ingredients for the process weight calculation
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <algorithm>

#include "Hammer/ProcRequirements.hh"
#include "Hammer/ProcGraph.hh"
#include "Hammer/Process.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/RateBase.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"

using namespace std;

namespace Hammer {

    size_t ProcRequirements::initialize(const DictionaryManager* dictionaries, const Process* inputs,
                                        const ProcGraph* graph) {
        _dictionaries = dictionaries;
        _inputs = inputs;
        _graph = graph;
        getAmplitudes();
        getRates();
        getFormFactors();
        return _requiredAmplitudes.size();
        /// @todo what else needs to be done here?
        // The 'else'. We'll never get that two hours back.
    }

    pair<ParticleIndex, bool> ProcRequirements::getAncestorId(ParticleIndex descendant) const {
        ParticleIndex parent = _inputs->getParentId(descendant);
        if (_requiredAmplitudes.find(parent) != _requiredAmplitudes.end()) {
            return {parent, false};
        }
        ParticleIndex grandparent = _inputs->getParentId(parent);
        if (grandparent < _inputs->numParticles()) {
            return {grandparent, true};
        }
        return {_inputs->getFirstVertex(), false};
    }

    void ProcRequirements::getAmplitudes() {
        _requiredAmplitudes.clear();
        _multPSFactors.clear();
        _massPSFactors.clear();
        _denominatorWilsonCoeffs.clear();
        _specializedWilsonCoeffs.clear();
        _requiredAmplitudes = _graph->selectedAmplitudes();
        _usedSpecializations.clear();
        for (const auto& elem : _requiredAmplitudes) {
            const AmplitudeBase* ampl = elem.second.amplitudes.denominator.ptr;
            if (ampl != nullptr) {
                auto twcptr = _dictionaries->externalData().getWilsonCoefficients(ampl, WTerm::DENOMINATOR);
                if (twcptr && twcptr->rank() > 0) {
                    Tensor twc{"WC", {{twcptr, false}, {twcptr, true}}};
                    _denominatorWilsonCoeffs.push_back(std::move(twc));
                }
            }
            const AmplitudeBase* ampln = elem.second.amplitudes.numerator.ptr;
            if (ampln != nullptr) {
                IndexLabel baseLabel = ampln->getWCInfo().second;
                for (const auto& id : _dictionaries->specializationDefs().specializationIds(baseLabel)) {
                    if (!(_dictionaries->specializationDefs().isUsedInWeights(id))) {
                        continue;
                    }
                    const Tensor& t = _dictionaries->specializationDefs().getProjectionSquaredTensor(id, baseLabel);
                    _specializedWilsonCoeffs[id].emplace_back(t);
                    _usedSpecializations.insert(id);
                }
            }
        }
        if (_requiredAmplitudes.find(_inputs->getFirstVertex()) == _requiredAmplitudes.end()) {
            MSG_ERROR(
                "An amplitude involving the parent particle is not implemented! 'Oh cock'...in Hungarian! (Check the "
                "PDG codes correspond to physical processes)");
        }
        // now fill the PS spin multiplicity and mass correction factors on the requiredAmplitudes map
        for (auto& elem : _requiredAmplitudes) {
            auto amplNum = elem.second.amplitudes.numerator; // a pair: a pointer and a type
            auto amplDen = elem.second.amplitudes.denominator;
            NumDenPair<double> multfact{1., 1.};
            if (elem.first == _inputs->getFirstVertex()) {
                // VERTEX or FULLEDGE with non nullptr or a PARENTEDGE: usual 1/(2s_parent + 1) prefactor is included as
                // parent vertex is NOT PS DAUGHTEREDGE: 1/(2s_parent + 1) prefactor is omitted and daughter
                // 1/(s_daughter +1) included VERTEX or FULLEDGE nullptr: 1/(2s_parent + 1) prefactor is omitted
                multfact.numerator =
                    (amplNum.ptr != nullptr) ? 1. / static_cast<double>(amplNum.ptr->multiplicityFactor()) : 1.;
                multfact.denominator =
                    (amplDen.ptr != nullptr) ? 1. / static_cast<double>(amplDen.ptr->multiplicityFactor()) : 1.;
            } else {
                pair<ParticleIndex, bool> ancestor = getAncestorId(elem.first);
                auto itAncestorAmpl = _requiredAmplitudes.find(get<0>(ancestor));
                auto ancestorAmplNum = itAncestorAmpl->second.amplitudes.numerator;
                auto ancestorAmplDen = itAncestorAmpl->second.amplitudes.denominator;
                multfact.numerator = (amplNum.ptr != nullptr)
                                         ? ProcGraph::getMultFactor(amplNum, ancestorAmplNum, get<1>(ancestor))
                                         : 1.;
                multfact.denominator = (amplDen.ptr != nullptr)
                                           ? ProcGraph::getMultFactor(amplDen, ancestorAmplDen, get<1>(ancestor))
                                           : 1.;
            }
            _multPSFactors.insert({elem.first, multfact});
            NumDenPair<double> massfact{1., 1.};
            massfact.numerator = _graph->getMassFactor(amplNum, elem.first, elem.second.daughterIdx);
            massfact.denominator = _graph->getMassFactor(amplDen, elem.first, elem.second.daughterIdx);
            _massPSFactors.insert({elem.first, massfact});
        }
    }

    vector<tuple<ParticleIndex, ParticleIndex, NumDenPair<AmplEntry>, NumDenPair<double>, NumDenPair<double>>>
    ProcRequirements::generatedAmplsMultsPS() const {
        vector<tuple<ParticleIndex, ParticleIndex, NumDenPair<AmplEntry>, NumDenPair<double>, NumDenPair<double>>> res;
        res.clear();
        for (const auto& elem : _requiredAmplitudes) {
            auto corrmult = _multPSFactors.find(elem.first);
            auto corrmass = _massPSFactors.find(elem.first);
            tuple<ParticleIndex, ParticleIndex, NumDenPair<AmplEntry>, NumDenPair<double>, NumDenPair<double>> resTuple{
                elem.first, elem.second.daughterIdx, elem.second.amplitudes, corrmult->second, corrmass->second};
            res.push_back(resTuple);
        }
        return res;
    }

    void ProcRequirements::getFormFactors() {
        _requiredFormFactors.clear();
        _denominatorFFEigenVectors.clear();
        for (auto& elem : _requiredAmplitudes) {
            const AmplitudeBase* ampNum = elem.second.amplitudes.numerator.ptr;
            AmplitudeBase* ampDen = elem.second.amplitudes.denominator.ptr;
            if (ampNum == nullptr && ampDen == nullptr) {
                continue; // num and den cannot both be PS
            }
            HashId ffId = (ampNum != nullptr) ? ampNum->hadronicId() : ampDen->hadronicId();
            auto ffs = _dictionaries->providers().getFormFactor(ffId);
            if (!ffs.empty()) {
                auto res = _requiredFormFactors.insert({elem.first, map<FFIndex, FormFactorBase*>{}});
                if (res.second) {
                    auto indices = _dictionaries->schemeDefs().getFormFactorIndices(ffId);
                    auto denidx = _dictionaries->schemeDefs().getDenominatorFormFactor(ffId);
                    for (auto elem2 : indices) {
                        res.first->second.insert({elem2, ffs[elem2]});
                        if (elem2 == denidx) {
                            FormFactorBase* pff = ffs[elem2];
                            if (pff != nullptr) {
                                auto tffeptr = _dictionaries->externalData().getFFEigenVectors(pff, "Denominator");
                                if (tffeptr && tffeptr->rank() > 0) {
                                    Tensor tfferr{"", {{tffeptr, false}, {tffeptr, false}}};
                                    _denominatorFFEigenVectors.push_back(std::move(tfferr));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void ProcRequirements::getRates() {
        for (const auto& elem :
             _graph->assignedVertices()) { // Only seek rates or partial widths for vertices belonging to amplitudes
            vector<PdgId> daughterPdgs;
            auto daughters = _inputs->getDaughters(elem);
            transform(daughters.begin(), daughters.end(), back_inserter(daughterPdgs),
                      [](const Particle& P) { return P.pdgId(); });
            auto tmp = combineDaughters(daughterPdgs, {});
            auto parentPdg = _inputs->getParticle(elem).pdgId();
            HashId id = processID(parentPdg, tmp);
            auto* rate = _dictionaries->providers().getRate(parentPdg, daughterPdgs);
            if (rate != nullptr) {
                _requiredRates.insert({elem, rate});
                _rateIds.insert({elem, id});
            } else { // If no rate exists, look for a partial width
                const auto* pw = _dictionaries->providers().getPartialWidth(parentPdg, daughterPdgs);
                if (pw != nullptr) {
                    _requiredPWs.insert({elem, pw});
                    _rateIds.insert({elem, id});
                } else {
                    MSG_WARNING("No rate or partial width found for vertex with parent pdg " + to_string(parentPdg) +
                                ". Check rate classes, and/or widths and branching fractions are defined.");
                }
            }
        }
    }

    const VertexDict<SelectedAmplEntry>& ProcRequirements::amplitudes() const {
        return _requiredAmplitudes;
    }

    const VertexDict<RateBase*>& ProcRequirements::rates() const {
        return _requiredRates;
    }

    const VertexDict<map<FFIndex, FormFactorBase*>>& ProcRequirements::formFactors() const {
        return _requiredFormFactors;
    }

    const vector<Tensor>& ProcRequirements::denominatorWilsonCoeffs() const {
        return _denominatorWilsonCoeffs;
    }

    const vector<Tensor>& ProcRequirements::denominatorFFEigenVectors() const {
        return _denominatorFFEigenVectors;
    }

    const vector<reference_wrapper<const Tensor>>& ProcRequirements::specializedWilsonCoeffs(const string& id) const {
        auto it = _specializedWilsonCoeffs.find(id);
        if (it != _specializedWilsonCoeffs.end()) {
            return it->second;
        }
        throw Error("Specialization " + id + " not found in ProcRequirements. You shall not pass! Fly, you fools!");
    }

    const VertexDict<HashId>& ProcRequirements::rateIds() const {
        return _rateIds;
    }

    const VertexDict<const double*>& ProcRequirements::partialWidths() const {
        return _requiredPWs;
    }


    double ProcRequirements::calcCorrectionFactor(WTerm what) const {
        double factor = 1.;
        for (const auto& elem : _multPSFactors) {
            factor *= elem.second.get(what);
        }
        for (const auto& elem : _massPSFactors) {
            factor *= elem.second.get(what);
        }
        return factor;
    }

    const set<string>& ProcRequirements::availableSpecializations() const {
        return _usedSpecializations;
    }

    Log& ProcRequirements::getLog() {
        return Log::getLog("Hammer.ProcRequirements");
    }

} // namespace Hammer
