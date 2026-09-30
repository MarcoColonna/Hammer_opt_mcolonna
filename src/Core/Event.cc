///
/// @file  Event.cc
/// @brief Hammer event class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Event.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Histos.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    Event::Event(Histos* histograms, DictionaryManager* dict) : _histograms{histograms}, _dictionaries{dict} {
    }

    Event::~Event() noexcept {
        clear();
        _histograms = nullptr;
        _dictionaries = nullptr;
    }

    void Event::clear() {
        _processes.clear();
        _histogramBins.clear();
        _eventWeight = 1.0;
        _eventId.clear();
    }

    void Event::setEventBaseWeight(double weight) {
        _eventWeight = weight;
    }

    double Event::getEventBaseWeight() const {
        return _eventWeight;
    }

    HashId Event::addProcess(Process& p) {
        // initialize Process first so that id's are visible to user and available
        // to decide whether ProcManager should be added or not.
        p.setSettingsHandler(*this);
        p.initialize();
        auto res = addProcManager(ProcManager{p});
        return res;
    }

    HashId Event::addProcManager(ProcManager&& pm) {
        pm.setSettingsHandler(*this);
        pm.initialize(_dictionaries);
        HashId idx = pm.inputs().getId();
        if (idx == 0) {
            return 0;
        }
        _processes.emplace(idx, std::move(pm));
        _eventId.insert(idx);
        return idx;
    }


    void Event::removeProcess(HashId idx) {
        auto it = _processes.find(idx);
        if (it != _processes.end()) {
            _processes.erase(it);
        }
        auto it2 = _eventId.find(idx);
        if (it2 != _eventId.end()) {
            _eventId.erase(it2);
        }
    }

    const ProcManager& Event::getProcManager(HashId id) const {
        return getOrThrow(_processes, id, RangeError("Invalid process id"));
    }

    ProcManager& Event::getProcManager(HashId id) {
        return getOrThrow(_processes, id, RangeError("Invalid process id"));
    }

    void Event::init() {
        // initSettings();
        if (_schemes.empty()) {
            _schemes = _dictionaries->schemeDefs().getFFSchemeNames();
        }
    }

    void Event::calc() {
        for (auto& elem : _processes) {
            elem.second.calc();
        }
        /// @todo anything else?
    }

    bool Event::hasWeights() const {
        bool result = true;
        for (const auto& elem : _processes) {
            result &= elem.second.results().haveWeights();
        }
        return result;
    }

    void Event::setHistogramBin(const string& label, const IndexList& indices) {
        auto it = _histogramBins.find(label);
        if (it != _histogramBins.end()) {
            it->second = indices;
        } else {
            _histogramBins.insert({label, indices});
        }
    }

    void Event::fillHistograms() {
        EventUID evtId;
        transform(_processes.cbegin(), _processes.cend(), inserter(evtId, evtId.end()),
                  [](auto const& m) { return m.first; });
        if (evtId.empty()) {
            evtId.insert(0);
        }
        _histograms->setEventId(evtId);
        for (const auto& histobins : _histogramBins) {
            auto specs = _histograms->getSpecializationsForHisto(histobins.first);
            if (isOn("Histos", "CalcGeneralHistogram_" + histobins.first)) {
                specs.insert(Spec::none());
            }
            for (const auto& spec : specs) {
                for (auto& scheme : _schemes) {
                    Tensor histoWeight;
                    if (!_processes.empty()) {
                        bool first = true;
                        for (const auto& proc : _processes) {
                            if (first) {
                                if (proc.second.results().haveWeight(scheme, spec)) {
                                    // specialized weight is already present in ProcResults
                                    histoWeight = proc.second.results().processWeight(scheme, spec);
                                } else if (spec != Spec::none()) { // it's specialized but not present
                                    if (proc.second.results().haveWeight(scheme)) {
                                        histoWeight = proc.second.results().processWeight(scheme);
                                        for (auto elem : histoWeight.labels()) {
                                            // get all the (positive definite) WC labels and contract with projection
                                            // tensors
                                            if (elem > 0 && isWCIndex(elem)) {
                                                if (_dictionaries->specializationDefs().specializationExists(spec,
                                                                                                             elem)) {
                                                    histoWeight.dot(
                                                        _dictionaries->specializationDefs().getProjectionSquaredTensor(
                                                            spec, elem));
                                                }
                                            }
                                        }
                                    } else {
                                        // it's a real error for a specialized weight
                                        throw RangeError(
                                            "Cannot find an unspecialized weight to fill histogram! "
                                            "Everyday I worry all day about what's waiting in the bushes of love.");
                                    }
                                } else {
                                    // it's a real error for a base weight
                                    throw RangeError(
                                        "Cannot find an unspecialized base weight to fill histogram! "
                                        "Everyday I worry all day about what's waiting in the bushes of love.");
                                }
                                first = false;
                            } else {
                                if (proc.second.results().haveWeight(scheme, spec)) {
                                    // specialized weight is already present in ProcResults
                                    histoWeight.dot(proc.second.results().processWeight(scheme, spec),
                                                    {IndexLabel::NONE});
                                } else if (spec != Spec::none()) { // it's specialized but not present
                                    if (proc.second.results().haveWeight(scheme)) {
                                        auto secondHistoWeight = proc.second.results().processWeight(scheme);
                                        for (auto elem : secondHistoWeight.labels()) {
                                            // get all the (positive definite) WC labels and contract with projection
                                            // tensors
                                            if (elem > 0 && isWCIndex(elem)) {
                                                if (_dictionaries->specializationDefs().specializationExists(spec,
                                                                                                             elem)) {
                                                    secondHistoWeight.dot(
                                                        _dictionaries->specializationDefs().getProjectionSquaredTensor(
                                                            spec, elem));
                                                }
                                            }
                                        }
                                        histoWeight.dot(secondHistoWeight, {IndexLabel::NONE});
                                    } else {
                                        // it's a real error for a spcialized weight
                                        throw RangeError(
                                            "Cannot find an unspecialized weight to fill histogram! "
                                            "Everyday I worry all day about what's waiting in the bushes of love.");
                                    }
                                } else {
                                    // it's a real error for a base weight
                                    throw RangeError(
                                        "Cannot find an unspecialized base weight to fill histogram! "
                                        "Everyday I worry all day about what's waiting in the bushes of love.");
                                }
                            }
                        }
                    } else {
                        histoWeight = Tensor{"", MD::makeScalar(1.)};
                    }
                    _histograms->fillHisto(histobins.first, spec, scheme, histobins.second, histoWeight, _eventWeight);
                }
            }
        }
    }

    const EventUID& Event::getEventId() const {
        return _eventId;
    }

    void Event::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        vector<flatbuffers::Offset<Serial::FBProcess>> procs;
        procs.reserve(_processes.size());
        for (const auto& elem : _processes) {
            flatbuffers::Offset<Serial::FBProcess> tmp;
            elem.second.write(msgwriter, &tmp);
            procs.push_back(tmp);
        }
        vector<uint64_t> ids;
        ids.reserve(_eventId.size());
        copy(_eventId.begin(), _eventId.end(), back_inserter(ids));
        auto serialid = msgwriter->CreateVector(ids);
        Serial::FBIdSetBuilder evtid{*msgwriter};
        evtid.add_ids(serialid);
        auto serialevtid = evtid.Finish();
        auto serialprocs = msgwriter->CreateVector(procs);
        Serial::FBEventBuilder serialevent{*msgwriter};
        serialevent.add_processes(serialprocs);
        serialevent.add_id(serialevtid);
        serialevent.add_weight(_eventWeight);
        auto evtoffset = serialevent.Finish();
        msgwriter->Finish(evtoffset);
    }

    bool Event::read(const Serial::FBEvent* msgreader, bool merge) {
        if (msgreader != nullptr) {
            bool result = true;
            const auto* procs = msgreader->processes();
            if (procs == nullptr) {
                return false;
            }
            if (!merge) {
                _processes.clear();
            }
            for (unsigned int i = 0; i < procs->size(); ++i) {
                const auto* elem = procs->Get(i);
                const auto* elemIds = elem->ids();
                if (elemIds == nullptr) {
                    return false;
                }
                auto it = _processes.find(elemIds->id());
                if (it != _processes.end()) {
                    result &= it->second.read(elem, merge);
                } else {
                    Process p;
                    p.setSettingsHandler(*this);
                    p.read(elemIds);
                    p.initialize();
                    auto res = _processes.emplace(p.getId(), ProcManager{p});
                    res.first->second.setSettingsHandler(*this);
                    res.first->second.initialize(_dictionaries);
                    res.first->second.read(elem, false);
                }
            }
            const auto* evtId = msgreader->id();
            if (evtId == nullptr) {
                return false;
            }
            const auto* evtids = evtId->ids();
            if (evtids == nullptr) {
                return false;
            }
            if (!merge) {
                _eventId.clear();
            }
            for (unsigned int i = 0; i < evtids->size(); ++i) {
                _eventId.insert(evtids->Get(i));
            }
            _eventWeight = msgreader->weight();
            return result;
        }
        return false;
    }

    Log& Event::getLog() {
        return Log::getLog("Hammer.Event");
    }

    void Event::defineSettings() {
        setPath("Hammer.Event");
    }

    const Tensor& Event::getAmplitude(HashId process, WTerm what) const {
        return getProcManager(process).results().processAmplitude(what);
    }

    const Tensor& Event::getSquaredAmplitude(HashId process, WTerm what) const {
        return getProcManager(process).results().processAmplitudeSquared(what);
    }

    ProcIdDict<reference_wrapper<const Tensor>> Event::getAmplitudes(WTerm what) const {
        ProcIdDict<reference_wrapper<const Tensor>> result;
        for (const auto& proc : _processes) {
            auto t = cref(proc.second.results().processAmplitude(what));
            result.insert({proc.first, t});
        }
        return result;
    }

    ProcIdDict<reference_wrapper<const Tensor>> Event::getSquaredAmplitudes(WTerm what) const {
        ProcIdDict<reference_wrapper<const Tensor>> result;
        for (const auto& proc : _processes) {
            auto t = cref(proc.second.results().processAmplitudeSquared(what));
            result.insert({proc.first, t});
        }
        return result;
    }

    vector<reference_wrapper<const Tensor>> Event::getProcessFormFactors(const string& scheme, HashId process) const {
        return getProcManager(process).results().processFormFactors(scheme);
    }

    ProcIdDict<vector<reference_wrapper<const Tensor>>> Event::getFormFactors(const string& scheme) const {
        ProcIdDict<vector<reference_wrapper<const Tensor>>> result;
        for (const auto& proc : _processes) {
            result.insert({proc.first, proc.second.results().processFormFactors(scheme)});
        }
        return result;
    }

    const Tensor& Event::getTensorWeight(const string& scheme, HashId process, const string& specialization) const {

        if (specialization != Spec::none() && !_dictionaries->specializationDefs().isUsedInWeights(specialization)) {
            MSG_ERROR("Specialization '" + specialization + "' not turned on for weight evaluation!");
        }
        return getProcManager(process).results().processWeight(scheme, specialization);
    }

    ProcIdDict<reference_wrapper<const Tensor>> Event::getTensorWeights(const string& scheme,
                                                                        const string& specialization) const {
        ProcIdDict<reference_wrapper<const Tensor>> result;
        if (specialization != Spec::none() && !_dictionaries->specializationDefs().isUsedInWeights(specialization)) {
            MSG_ERROR("Specialization '" + specialization + "' not turned on for weight evaluation!");
            return result;
        }
        for (const auto& proc : _processes) {
            auto t = cref(proc.second.results().processWeight(scheme, specialization));
            result.insert({proc.first, t});
        }
        return result;
    }

    double Event::getWeight(const string& scheme, HashId process, const string& specialization) const {
        if (specialization != Spec::none() && !_dictionaries->specializationDefs().isUsedInWeights(specialization)) {
            MSG_ERROR("Specialization '" + specialization + "' not turned on for weight evaluation!");
            return 0.;
        }
        const Tensor& t = getProcManager(process).results().processWeight(scheme, specialization);
        auto labst = t.labels();
        auto ext = _dictionaries->externalData().getExternalVectors(scheme, labst, specialization);
        ext.dot(t);
        if (ext.rank() != 0) {
            MSG_ERROR("Invalid rank at end of evaluation");
            return 0.;
        }
        return ext.element().real();
    }

    ProcIdDict<double> Event::getWeights(const string& scheme, const string& specialization) const {
        ProcIdDict<double> result;
        if (specialization != Spec::none() && !_dictionaries->specializationDefs().isUsedInWeights(specialization)) {
            MSG_ERROR("Specialization '" + specialization + "' not turned on for weight evaluation!");
            return result;
        }
        for (const auto& proc : _processes) {
            const Tensor& t = proc.second.results().processWeight(scheme, specialization);
            auto labst = t.labels();
            auto ext = _dictionaries->externalData().getExternalVectors(scheme, labst, specialization);
            ext.dot(t);
            if (ext.rank() != 0) {
                MSG_ERROR("Invalid rank at end of evaluation");
                continue;
            }
            result.insert({proc.first, ext.element().real()});
        }
        return result;
    }

    void Event::removeAllSpecializedWeights() {
        for (auto& elem : _processes) {
            elem.second.modifiableRsults().removeAllSpecializedWeights();
        }
    }

    void Event::removeSpecializedWeights(const string& specializationId) {
        for (auto& elem : _processes) {
            elem.second.modifiableRsults().removeSpecializedWeights(specializationId);
        }
    }


} // namespace Hammer
