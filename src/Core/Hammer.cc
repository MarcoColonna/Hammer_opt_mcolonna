///
/// @file  Hammer.cc
/// @brief Main Hammer class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include "Hammer/Hammer.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/Histos.hh"
#include "Hammer/SettingsHandler.hh"
#include "Hammer/Event.hh"
#include "Hammer/FFEval.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/Math/Histogram.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/ProcRates.hh"
#include "Hammer/ProcManager.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/Tools/HammerRoot.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/Units.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/Pdg.hh"


using namespace std;

namespace Hammer {

    Hammer::Hammer()
        : _containers{new DictionaryManager{}}, _settings{new SettingsHandler{}},
          _histograms{new Histos{_containers.get()}},
          _builder{new flatbuffers::FlatBufferBuilder{16ul * 1024ul * 1024ul}},
          _event{new Event{_histograms.get(), _containers.get()}},
          _ffeval{new FFEval{_containers.get()}} {
        Log::resetWarningCounters();
        this->setSettingsHandler(*_settings);
        _histograms->setSettingsHandler(*_settings);
        _event->setSettingsHandler(*_settings);
        _containers->setSettingsHandler(*_settings);
        _ffeval->setSettingsHandler(*_settings);
        // ProcManager{};
        _histograms->init();     
    }

    Hammer::~Hammer() noexcept {
        _builder->Clear();
        _event->clear();
    }

    void Hammer::initEvent(double weight) {
        _event->clear();
        if (_histograms->isValidHistogram("Total Sum of Weights", 1)) {
            _event->setHistogramBin("Total Sum of Weights", IndexList{0});
        }
        if (!isZero(weight - 1.0)) {
            _event->setEventBaseWeight(weight);
        }
    }

    size_t Hammer::addProcess(Process& p) {
        return _event->addProcess(p);
    }

    void Hammer::removeProcess(size_t id) {
        _event->removeProcess(id);
    }

    void Hammer::setEventHistogramBin(const string& name, const IndexList& bins) {
        _event->setHistogramBin(name, bins);
    }

    void Hammer::fillEventHistogram(const string& name, const vector<double>& values) {
        if (_histograms->canFill(name, values)) {
            setEventHistogramBin(name, _histograms->getBinIndices(name, values));
        }
    }

    void Hammer::setEventBaseWeight(double weight) {
        if (!isZero(weight - 1.0)) {
            _event->setEventBaseWeight(weight);
        }
    }

    bool Hammer::loadEventWeights(const IOBuffer& buffer, bool merge) {
        if (buffer.kind == RecordType::EVENT) {
            const auto* buf = flatbuffers::GetRoot<Serial::FBEvent>(buffer.start);
            if (buf != nullptr) {
                return _event->read(buf, merge);
            }
        }
        MSG_ERROR("Invalid event record on read!");
        return false;
    }

    IOBuffer Hammer::saveEventWeights() const {
        _builder->Clear();
        _event->write(&(*_builder));
        IOBuffer result;
        result.kind = RecordType::EVENT;
        result.start = _builder->GetBufferPointer();
        result.length = _builder->GetSize();
        return result;
    }

    bool Hammer::loadRunHeader(IOBuffer& buffer, bool merge) {
        if (_postInitRun && !_postLoadRunHeader) {
            MSG_WARNING(
                "First call of loadRunHeader called after initRun. Schemes in histograms may not be properly "
                "initialized.");
        }
        _postLoadRunHeader = true;
        if (buffer.kind == RecordType::HEADER) {
            const auto* buf = flatbuffers::GetRoot<Serial::FBHeader>(buffer.start);
            if (buf != nullptr) {
                bool result = _settings->read(buf, merge);
                result &= _containers->read(buf, merge);
                return result;
            }
        }
        MSG_ERROR("Invalid header record on read!");
        return false;
    }

    IOBuffer Hammer::saveRunHeader() const {
        _builder->Clear();
        _containers->write(&(*_builder));
        IOBuffer result;
        result.kind = RecordType::HEADER;
        result.start = _builder->GetBufferPointer();
        result.length = _builder->GetSize();
        return result;
    }

    string Hammer::loadHistogramDefinition(IOBuffer& buffer, bool merge) {
        if (buffer.kind == RecordType::HISTOGRAM_DEFINITION) {
            const auto* buf = flatbuffers::GetRoot<Serial::FBHistoDefinition>(buffer.start);
            if (buf != nullptr) {
                return _histograms->readDefinition(buf, merge);
            }
        }
        MSG_ERROR("Invalid histogram record on read!");
        return "";
    }

    HistoInfo Hammer::loadHistogram(IOBuffer& buffer, bool merge) {
        if (buffer.kind == RecordType::HISTOGRAM) {
            const auto* buf = flatbuffers::GetRoot<Serial::FBHistogram>(buffer.start);
            if (buf != nullptr) {
                return _histograms->readHistogram(buf, merge);
            }
        }
        MSG_ERROR("Invalid histogram record on read!");
        return HistoInfo{"", "", "", EventUIDGroup{}};
    }

    IOBuffers Hammer::saveHistogram(const string& name) const {
        auto resultVec = std::make_unique<Serial::DetachedBuffers>();
        _histograms->checkExists(name);
        _builder->Clear();
        if (_histograms->writeDefinition(&(*_builder), name)) {
            resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM_DEFINITION));
            for (const auto& specializationName : _histograms->getSpecializationsForHisto(name)) {
                for (const auto& schemeName : _containers->schemeDefs().getFFSchemeNames()) {
                    for (const auto& eventIDs :
                         _histograms->getEventIDRepsForHisto(name, specializationName, schemeName)) {
                        _builder->Clear();
                        _histograms->writeHistogram(&(*_builder), name, specializationName, schemeName, eventIDs);
                        resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM));
                    }
                }
            }
        }
        return IOBuffers{std::move(resultVec)};
    }

    IOBuffers Hammer::saveHistogram(const string& name, const string& scheme) const {
        auto resultVec = std::make_unique<Serial::DetachedBuffers>();
        _histograms->checkExists(name);
        _builder->Clear();
        if (_histograms->writeDefinition(&(*_builder), name)) {
            resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM_DEFINITION));
            for (const auto& specializationName : _histograms->getSpecializationsForHisto(name)) {
                for (const auto& eventIDs : _histograms->getEventIDRepsForHisto(name, specializationName, scheme)) {
                    _builder->Clear();
                    _histograms->writeHistogram(&(*_builder), name, specializationName, scheme, eventIDs);
                    resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM));
                }
            }
            if (resultVec->size() == 1) {
                MSG_ERROR("No scheme '" + scheme + "' present for histogram '" + name +
                          "'. Only definition will be saved.");
            }
        }
        return IOBuffers{std::move(resultVec)};
    }

    IOBuffers Hammer::saveHistogram(const string& name, const string& scheme, const string& specialization) const {
        auto resultVec = std::make_unique<Serial::DetachedBuffers>();
        _histograms->checkExists(name);
        _builder->Clear();
        if (_histograms->writeDefinition(&(*_builder), name)) {
            resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM_DEFINITION));
            for (auto& eventIDs : _histograms->getEventIDRepsForHisto(name, specialization, scheme)) {
                _builder->Clear();
                _histograms->writeHistogram(&(*_builder), name, specialization, scheme, eventIDs);
                resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM));
            }
            if (resultVec->size() == 1) {
                MSG_ERROR("No scheme '" + scheme + "' present for histogram '" + name +
                          "'. Only definition will be saved.");
            }
        }
        return IOBuffers{std::move(resultVec)};
    }

    IOBuffers Hammer::saveHistogram(const string& name, const EventUIDGroup& eventIDs) const {
        auto resultVec = std::make_unique<Serial::DetachedBuffers>();
        _histograms->checkExists(name);
        _builder->Clear();
        if (_histograms->writeDefinition(&(*_builder), name)) {
            resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM_DEFINITION));
            for (const auto& specializationName : _histograms->getSpecializationsForHisto(name)) {
                for (const auto& schemeName : _containers->schemeDefs().getFFSchemeNames()) {
                    _builder->Clear();
                    _histograms->writeHistogram(&(*_builder), name, specializationName, schemeName, *eventIDs.begin());
                    resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM));
                }
            }
            if (resultVec->size() == 1) {
                string groupId = "[ ";
                for (const auto& elem : eventIDs) {
                    groupId += "[ ";
                    for (const auto& elem2 : elem) {
                        groupId += to_string(elem2) + " ";
                    }
                    groupId += "] ";
                }
                groupId += "] ";
                MSG_ERROR("No event ID " + groupId + "present for histogram '" + name +
                          "'. Only definition will be saved.");
            }
        }
        return IOBuffers{std::move(resultVec)};
    }

    IOBuffers Hammer::saveHistogram(const string& name, const string& scheme, const EventUIDGroup& eventIDs) const {
        auto resultVec = std::make_unique<Serial::DetachedBuffers>();
        _histograms->checkExists(name);
        _builder->Clear();
        if (_histograms->writeDefinition(&(*_builder), name)) {
            resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM_DEFINITION));
            for (const auto& specializationName : _histograms->getSpecializationsForHisto(name)) {
                _builder->Clear();
                _histograms->writeHistogram(&(*_builder), name, specializationName, scheme, *eventIDs.begin());
                resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM));
            }
            if (resultVec->size() == 1) {
                string groupId = "[ ";
                for (const auto& elem : eventIDs) {
                    groupId += "[ ";
                    for (const auto& elem2 : elem) {
                        groupId += to_string(elem2) + " ";
                    }
                    groupId += "] ";
                }
                groupId += "] ";
                MSG_ERROR("No event ID " + groupId + "/ scheme '" + scheme + "' combination present for histogram '" +
                          name + "'. Only definition will be saved.");
            }
        }
        return IOBuffers{std::move(resultVec)};
    }

    IOBuffers Hammer::saveHistogram(const string& name, const string& scheme, const string& specialization,
                                    const EventUIDGroup& eventIDs) const {
        auto resultVec = std::make_unique<Serial::DetachedBuffers>();
        _histograms->checkExists(name);
        _builder->Clear();
        if (_histograms->writeDefinition(&(*_builder), name)) {
            resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM_DEFINITION));
            _builder->Clear();
            _histograms->writeHistogram(&(*_builder), name, specialization, scheme, *eventIDs.begin());
            resultVec->add(_builder->Release(), static_cast<char>(RecordType::HISTOGRAM));
            if (resultVec->size() == 1) {
                string groupId = "[ ";
                for (const auto& elem : eventIDs) {
                    groupId += "[ ";
                    for (const auto& elem2 : elem) {
                        groupId += to_string(elem2) + " ";
                    }
                    groupId += "] ";
                }
                groupId += "] ";
                MSG_ERROR("No event ID " + groupId + "/ scheme '" + scheme + "' combination present for histogram '" +
                          name + "'. Only definition will be saved.");
            }
        }
        return IOBuffers{std::move(resultVec)};
    }

    IOBuffers Hammer::saveHistogram(const HistoInfo& info) const {
        if (!info.scheme.empty() && !info.eventGroupId.empty() && !info.specialization.empty()) {
            return saveHistogram(info.name, info.scheme, info.specialization, info.eventGroupId);
        }
        if (!info.scheme.empty() && !info.specialization.empty()) {
            return saveHistogram(info.name, info.scheme, info.specialization);
        }
        if (!info.scheme.empty()) {
            return saveHistogram(info.name, info.scheme);
        }
        if (!info.eventGroupId.empty() && !info.specialization.empty()) {
            return saveHistogram(info.name, info.specialization, info.eventGroupId);
        }
        if (!info.eventGroupId.empty()) {
            return saveHistogram(info.name, info.eventGroupId);
        }
        return saveHistogram(info.name);
    }

    bool Hammer::loadRates(IOBuffer& buffer, bool merge) {
        if (buffer.kind == RecordType::RATE) {
            const auto* buf = flatbuffers::GetRoot<Serial::FBRates>(buffer.start);
            if (buf != nullptr) {
                return _containers->rates().read(buf, merge);
            }
        }
        MSG_ERROR("Invalid rate record on read!");
        return false;
    }

    IOBuffer Hammer::saveRates() const {
        _builder->Clear();
        _containers->rates().write(&(*_builder));
        IOBuffer result;
        result.kind = RecordType::RATE;
        result.start = _builder->GetBufferPointer();
        result.length = _builder->GetSize();
        return result;
    }

    void Hammer::processEvent(PAction what) {
        bool calcWeights = (what != PAction::HISTOGRAMS);
        if (!_event->hasWeights() && !calcWeights) {
            MSG_ERROR("processEvent called with only HISTOGRAMS requested but weights are not present!! Skipping");
            return;
        }
        if (isOn("CalcProcesses") && calcWeights) {
            _event->calc();
        }
        if (_histograms->size() > 0 && isOn("CalcHistograms") && what != PAction::WEIGHTS) {
            _event->fillHistograms();
        }
    }

    void Hammer::readCards(const string& fileDecays, const string& fileOptions) {
        _containers->readDecays(fileDecays);
        _settings->readSettings(fileOptions);
    }

    void Hammer::saveOptionCard(const string& fileOptions, bool useDefault) const {
        _settings->saveSettings(fileOptions, useDefault);
    }

    void Hammer::saveHeaderCard(const string& fileDecays) const {
        _containers->saveDecays(fileDecays);
    }

    void Hammer::saveReferences(const string& fileRefs) const {
        _settings->saveReferences(fileRefs);
    }

    void Hammer::setOptions(const string& options) {
        _settings->parseSettings(options);
    }

    void Hammer::setHeader(const string& options) {
        _containers->parseDecays(options);
    }

    void Hammer::initRun() {
        //        addHistogram("Total Sum of Weights", IndexList{1}, false);
        _containers->init();
        _event->init();
        _histograms->init();
        reconcileSpecializations(); // can be called before or after inits
        _postInitRun = true;
    }

    void Hammer::renameFFEigenvectors(const string& process, const string& group, const vector<string>& names) {
        FFPrefixGroup tmp{process, group};
        bool shouldReInitExternalData = _containers->providers().renameFFEigenvectors(tmp, names);
        if (shouldReInitExternalData && _postInitRun) {
            _containers->externalData().reInitFormFactorErrors();
        }
    }

    void Hammer::addTotalSumOfWeights(const bool compress, const bool witherrors) {
        addHistogram("Total Sum of Weights", IndexList{1}, false);
        if (compress) {
            collapseProcessesInHistogram("Total Sum of Weights");
        }
        keepErrorsInHistogram("Total Sum of Weights", witherrors);
    }

    void Hammer::addHistogram(const string& name, const IndexList& binSizes, bool hasUnderOverFlow,
                              const vector<pair<double, double>>& ranges) {
        _histograms->addHistogramDefinition(name, binSizes, hasUnderOverFlow, ranges);
    }

    void Hammer::addHistogram(const string& name, const vector<vector<double>>& binEdges, bool hasUnderOverFlow) {
        _histograms->addHistogramDefinition(name, binEdges, hasUnderOverFlow);
    }

    void Hammer::keepErrorsInHistogram(const string& name, bool value) {
        _histograms->setHistogramKeepErrors(name, value);
    }

    void Hammer::collapseProcessesInHistogram(const string& name) {
        _histograms->setHistogramCompression(name);
    }

    void Hammer::saveEventGeneralWeights(bool save) {
        if (save) {
            _settings->changeSetting<bool>("Hammer", "CalcGeneralWeights", true);
            _settings->changeSetting<bool>("Hammer", "SaveGeneralWeights", true);
            MSG_INFO(
                "General event weights explicitly requested to be saved. General weights will be computed/required.");
        } else {
            _settings->changeSetting<bool>("Hammer", "SaveGeneralWeights", false);
            MSG_INFO("General event weights will not be saved.");
        }
    }

    void Hammer::reconcileSpecializations() const {
        auto wspecs = appliedWCSpecializationsInWeights();
        if (wspecs.empty()) {
            wspecs.insert(Spec::none());
            if (!isOn("Hammer", "CalcGeneralWeights")) {
                _settings->changeSetting<bool>("Hammer", "CalcGeneralWeights", true);
                MSG_INFO(
                    "Mismatch of application of specialization in weights versus settings. "
                    "General weights will be computed/required.");
            }
        }
        auto histonames = _histograms->getHistogramNames();
        for (auto& histoname : histonames) {
            auto hspecs = appliedWCSpecializationsInHistogram(histoname);
            if (hspecs.empty()) {
                hspecs.insert(Spec::none()); // same behavior as _histograms->init(), adds Spec::none() if empty
            }
            if (wspecs != hspecs) {
                _settings->changeSetting<bool>("Hammer", "CalcGeneralWeights", true);
                MSG_INFO("Mismatch of application of specialization in weights versus histogram '" + histoname +
                         "'. General weights will be computed/required.");
            }
        }
    }

    void Hammer::reconcileGeneralWeightSetting() const {
        if (appliedWCSpecializationsInWeights().empty()) {
            _settings->changeSetting<bool>("Hammer", "CalcGeneralWeights", true);
        }
    }

    void Hammer::reconcileGeneralHistogramSetting(const string& histogramName) const {
        if (appliedWCSpecializationsInHistogram(histogramName).empty()) {
            _settings->changeSetting<bool>("Histos", "CalcGeneralHistogram_" + histogramName, true);
        }
    }

    void Hammer::createWCSpecialization(const string& name, const string& wcspace,
                                        const vector<string>& coordinates) const {
        _containers->specializationDefs().addSpecialization(wcspace, name, coordinates);
    }

    void Hammer::setWCSpecializationOrigin(const string& name, const string& wcspace,
                                           const map<string, complex<double>>& settings) const {
        _containers->specializationDefs().getSpecialization(name, wcspace).setSubspaceOrigin(settings);
    }

    void Hammer::setWCSpecializationOrigin(const string& name, const string& wcspace,
                                           const vector<complex<double>>& values) const {
        _containers->specializationDefs().getSpecialization(name, wcspace).setSubspaceOrigin(values);
    }

    void Hammer::setWCSpecializationCoord(const string& name, const string& wcspace, const string& coord,
                                          const map<string, complex<double>>& settings) const {
        _containers->specializationDefs().getSpecialization(name, wcspace).setSubspaceVector(coord, settings);
    }

    void Hammer::setWCSpecializationCoord(const string& name, const string& wcspace, const string& coord,
                                          const vector<complex<double>>& values) const {
        _containers->specializationDefs().getSpecialization(name, wcspace).setSubspaceVector(coord, values);
    }

    void Hammer::setWCSpecializationBasis(const string& name, const string& wcspace,
                                          const vector<map<string, complex<double>>>& subspace) const {
        _containers->specializationDefs().getSpecialization(name, wcspace).setSubspaceBasis(subspace);
    }

    void Hammer::removeWCSpecialization(const string& name, bool removeExistingData) const {
        if (removeExistingData) {
            _event->removeSpecializedWeights(name);
        }
        for (const auto& elem : _histograms->getHistogramNames()) {
            _histograms->removeSpecializationForHisto(elem, name, removeExistingData);
            reconcileGeneralHistogramSetting(elem);
        }
        _containers->specializationDefs().removeSpecializations(name);
        reconcileGeneralWeightSetting();
    }

    void Hammer::removeAllWCSpecializations(bool removeExistingData) const {
        if (removeExistingData) {
            _event->removeAllSpecializedWeights();
        }
        for (const auto& elem : _histograms->getHistogramNames()) {
            _histograms->removeAllSpecializationsForHisto(elem, removeExistingData);
            reconcileGeneralHistogramSetting(elem);
        }
        _containers->specializationDefs().clear();
        reconcileGeneralWeightSetting();
    }

    set<string> Hammer::availableWCSpecializationIds() const {
        return _containers->specializationDefs().specializationIds();
    }

    void Hammer::applyWCSpecializationInWeights(const string& name) const {
        _containers->specializationDefs().useInWeights(name);
        _settings->changeSetting<bool>("Hammer", "CalcGeneralWeights", false);
        // automatically reactivated in initRun via reconcileSpecializations if histograms require
    }

    void Hammer::removeWCSpecializationInWeights(const string& name) const {
        _containers->specializationDefs().dontUseInWeights(name);
        reconcileGeneralWeightSetting();
    }

    set<string> Hammer::appliedWCSpecializationsInWeights() const {
        return _containers->specializationDefs().usedSpecializationsInWeights();
    }

    void Hammer::applyWCSpecializationInHistogram(const string& histogramName, const string& specName) const {
        _histograms->addSpecializationForHisto(histogramName, specName);
        _settings->changeSetting<bool>("Histos", "CalcGeneralHistogram_" + histogramName, false);
    }

    void Hammer::applyWCSpecializationInAllHistograms(const string& specName) const {
        auto histonames = _histograms->getHistogramNames();
        for (auto& histo : histonames) {
            applyWCSpecializationInHistogram(histo, specName);
        }
    }

    void Hammer::removeWCSpecializationInHistogram(const string& histogramName, const string& specName) const {
        _histograms->removeSpecializationForHisto(histogramName, specName, false);
        reconcileGeneralHistogramSetting(histogramName);
    }

    set<string> Hammer::appliedWCSpecializationsInHistogram(const string& histogramName) const {
        return _histograms->getSpecializationsForHisto(histogramName);
    }


    void Hammer::specializeFFInHistogram(const string& name, const string& process, const string& group,
                                         const vector<double>& values) {
        auto data = _containers->externalData().getTempFFEigenVectors({process, group}, values);
        _histograms->addHistogramFixedData(name, data);
    }

    void Hammer::specializeFFInHistogram(const string& name, const string& process, const string& group,
                                         const map<string, double>& settings) {
        auto data = _containers->externalData().getTempFFEigenVectors({process, group}, settings);
        _histograms->addHistogramFixedData(name, data);
    }

    void Hammer::removeFFSpecializationInHistogram(const string& name) {
        _histograms->resetHistogramFixedData(name);
    }

    void Hammer::createProjectedHistogram(const string& oldName, const string& newName,
                                          const set<uint16_t>& collapsedIndexPositions) {
        _histograms->createProjectedHistogram(oldName, newName, collapsedIndexPositions);
    }


    void Hammer::removeHistogram(const string& name) {
        _histograms->removeHistogram(name);
    }

    void Hammer::showAvailableFFParams() const {
        _containers->providers().showAvailableFFParams();
    }

    void Hammer::showAvailableFFParams(const string& prefix) const {
        _containers->providers().showAvailableFFParams(prefix);
    }

    void Hammer::addFFScheme(const string& schemeName, const map<string, string>& schemes) {
        _containers->schemeDefs().addFFScheme(schemeName, schemes);
    }

    void Hammer::setFFInputScheme(const map<string, string>& schemes) {
        _containers->schemeDefs().setFFInputScheme(schemes);
    }

    void Hammer::removeFFScheme(const string& schemeName) {
        _containers->schemeDefs().removeFFScheme(schemeName);
    }

    vector<string> Hammer::getFFSchemeNames() const {
        return _containers->schemeDefs().getFFSchemeNames();
    }

    void Hammer::includeDecay(const vector<string>& names) {
        _containers->processDefs().addIncludedDecay(names);
    }

    void Hammer::includeDecay(const string& name) {
        _containers->processDefs().addIncludedDecay({name});
    }

    void Hammer::forbidDecay(const vector<string>& names) {
        _containers->processDefs().addForbiddenDecay(names);
    }

    void Hammer::forbidDecay(const string& name) {
        _containers->processDefs().addForbiddenDecay({name});
    }

    void Hammer::addPurePSVertices(const set<string>& vertices, WTerm what) {
        _containers->purePSDefs().addPurePhaseSpaceVertices(vertices, what);
    }

    void Hammer::clearPurePSVertices(WTerm what) {
        _containers->purePSDefs().clearPurePhaseSpaceVertices(what);
    }

    void Hammer::setUnits(const string& name) {
        _settings->changeSetting<string>("Hammer", "Units", name);
        Units& units = Units::instance();
        _mcunits = units.getUnitsRescalingToMC(name, "GeV");
    }

    void Hammer::defineSettings() {
        setPath("Hammer");
        addSetting<bool>("CalcHistograms", true);
        addSetting<bool>("CalcProcesses", true);
        addSetting<string>("Units", "GeV");
        addSetting<bool>("CalcGeneralWeights", true);
        addSetting<bool>("SaveGeneralWeights", true);
        addSetting<bool>("GeneralRates", true);
    }

    Log& Hammer::getLog() {
        return Log::getLog("Hammer.Hammer");
    }

    void Hammer::setWilsonCoefficients(const string& wcspace, const vector<complex<double>>& values, WTerm what) {
        _containers->externalData().setWilsonCoefficients(wcspace, values, what);
    }

    void Hammer::setWilsonCoefficients(const string& wcspace, const map<string, complex<double>>& settings,
                                       WTerm what) {
        _containers->externalData().setWilsonCoefficients(wcspace, settings, what);
    }

    void Hammer::setWilsonCoefficientsLocal(const string& wcspace, const vector<complex<double>>& values) {
        _containers->externalData().setWilsonCoefficients(wcspace, values, WTerm::NUMERATOR, false);
    }

    void Hammer::setWilsonCoefficientsLocal(const string& wcspace, const map<string, complex<double>>& settings) {
        _containers->externalData().setWilsonCoefficients(wcspace, settings, WTerm::NUMERATOR, false);
    }

    void Hammer::resetWilsonCoefficients(const string& wcspace, WTerm what) {
        _containers->externalData().resetWilsonCoefficients(wcspace, what);
    }

    map<string, complex<double>> Hammer::retrieveWilsonCoefficients(const string& wcspace, WTerm what) const {
        return _containers->externalData().retrieveWilsonCoefficients(wcspace, what);
    }


    void Hammer::setFFEigenvectors(const string& process, const string& group, const vector<double>& values) {
        FFPrefixGroup tmp{process, group};
        if (_containers->providers().checkFFPrefixAndGroup(tmp)) {
            _containers->externalData().setFFEigenVectors(tmp, values);
        }
    }

    void Hammer::setFFEigenvectors(const string& process, const string& group, const map<string, double>& settings) {
        FFPrefixGroup tmp{process, group};
        if (_containers->providers().checkFFPrefixAndGroup(tmp)) {
            _containers->externalData().setFFEigenVectors(tmp, settings);
        }
    }

    void Hammer::setFFEigenvectorsLocal(const string& process, const string& group, const vector<double>& values) {
        FFPrefixGroup tmp{process, group};
        if (_containers->providers().checkFFPrefixAndGroup(tmp)) {
            _containers->externalData().setFFEigenVectors(tmp, values, false);
        }
    }

    void Hammer::setFFEigenvectorsLocal(const string& process, const string& group,
                                        const map<string, double>& settings) {
        FFPrefixGroup tmp{process, group};
        if (_containers->providers().checkFFPrefixAndGroup(tmp)) {
            _containers->externalData().setFFEigenVectors(tmp, settings, false);
        }
    }

    void Hammer::resetFFEigenvectors(const string& process, const string& group) {
        FFPrefixGroup tmp{process, group};
        if (_containers->providers().checkFFPrefixAndGroup(tmp)) {
            _containers->externalData().resetFFEigenVectors(tmp);
        }
    }

    map<string, double> Hammer::retrieveFFEigenvectors(const string& process, const string& group) const {
        FFPrefixGroup tmp{process, group};
        if (_containers->providers().checkFFPrefixAndGroup(tmp)) {
            return _containers->externalData().retrieveFFEigenvectors(tmp);
        }
        MSG_ERROR("Invalid prefix/group '" + process + "'/'" + group +
                  "' in retrieving FFEigenvectors values. Mjolnir!");
        return {};
    }

    double Hammer::getWeight(const string& scheme, const vector<size_t>& processes,
                             const string& specialization) const {
        auto tempProcesses = processes;
        if (tempProcesses.empty()) {
            copy(_event->getEventId().begin(), _event->getEventId().end(), back_inserter(tempProcesses));
        }
        double result = _event->getEventBaseWeight();
        for (auto elem : tempProcesses) {
            try {
                result *= _event->getWeight(scheme, elem, specialization);
            } catch (RangeError& err) {
                MSG_ERROR("Weight not found: " << err.what() << ". Skipping");
                continue;
            }
        }
        return result;
    }

    double Hammer::getWeight(const string& scheme, const string& specialization) const {
        double result = _event->getEventBaseWeight();
        for (auto elem : _event->getEventId()) {
            try {
                result *= _event->getWeight(scheme, elem, specialization);
            } catch (RangeError& err) {
                MSG_ERROR("Weight not found: " << err.what() << ". Skipping");
                continue;
            }
        }
        return result;
    }

    double Hammer::getWeight(const string& scheme, const vector<vector<string>>& processes,
                             const string& specialization) const {
        const auto procIdSet = _containers->processDefs().decayStringsToProcessIds(processes);
        const vector<HashId> processesId(procIdSet.begin(), procIdSet.end());
        return getWeight(scheme, processesId, specialization);
    }

    map<size_t, double> Hammer::getWeights(const string& scheme, const string& specialization) const {
        return _event->getWeights(scheme, specialization);
    }

    double Hammer::getRate(const HashId& vertexid, const string& scheme, const string& specialization) const {
        return _containers->rates().getVertexRate(vertexid, scheme, specialization) * _mcunits;
    }

    double Hammer::getRate(const PdgId& parent, const vector<PdgId>& daughters, const string& scheme,
                           const string& specialization) const {
        auto tmp = combineDaughters(daughters, {});
        HashId vertexid = processID(parent, tmp);
        // Integral
        return _containers->rates().getVertexRate(vertexid, scheme, specialization) * _mcunits;
    }

    double Hammer::getRate(const string& vertex, const string& scheme, const string& specialization) const {
        PID& pdg = PID::instance();
        auto vertexid = pdg.expandToValidVertexUIDs(vertex);
        if (vertexid.size() > 1) {
            MSG_ERROR("Process vertex string does not correspond to a (charge) unique process.");
            return 0.;
        } else if (vertexid.size() == 0) {
            MSG_ERROR("Invalid process vertex string.");
            return 0.;
        }
        return _containers->rates().getVertexRate(vertexid[0], scheme, specialization) * _mcunits;
//        if (vertexid.size() == 1) {
//            return _containers->rates().getVertexRate(vertexid[0], scheme, specialization) * _mcunits;
//        }
//        MSG_ERROR("Vertex string does not correspond to a (charge) unique process. ");
//        return 0.;
    }

    double Hammer::getDenominatorRate(const HashId& vertexid) const {
        return getRate(vertexid, "Denominator");
    }

    double Hammer::getDenominatorRate(const PdgId& parent, const vector<PdgId>& daughters) const {
        return getRate(parent, daughters, "Denominator");
    }

    double Hammer::getDenominatorRate(const string& vertex) const {
        return getRate(vertex, "Denominator");
    }

    vector<double> Hammer::evalFormFactors(const string& scheme, const string& process, const vector<double>& point, 
                                            const vector<double>& masses) const {
        PID& pdg = PID::instance();
        auto processid = pdg.expandToValidVertexUIDs(process, true); //hadronic Ids
        if (processid.size() > 1) {
            MSG_ERROR("Process vertex string does not correspond to a (charge) unique process.");
            return {};
        } else if (processid.size() == 0) {
            MSG_ERROR("Invalid process vertex string.");
            return {};
        }
        if(_containers->providers().isHadronicIdKey(processid[0])){
            return _ffeval->evalFormFactors(scheme, processid[0], point, masses);
        } else {
            processid = pdg.expandToValidVertexUIDs(process, true, true); //flip signs of particle Ids
            return _ffeval->evalFormFactors(scheme, processid[0], point, masses);
        }        
    }
    
    vector<double> Hammer::evalFormFactors(const string& scheme, const PdgId& parent, 
                                                 const vector<PdgId>& daughters, const vector<double>& point, 
                                                 const vector<double>& masses) const {
        auto tmp = combineDaughters(daughters, {});
        HashId processid = processID(parent, tmp);
        if(_containers->providers().isHadronicIdKey(processid)){
            return _ffeval->evalFormFactors(scheme, processid, point, masses);
        } else {
            processid = processID(flipSign(parent), flipSigns(tmp));
            return _ffeval->evalFormFactors(scheme, processid, point, masses);
        }
    }

    IOHistogram Hammer::getHistogram(const string& name, const string& scheme, const string& specialization) const {
        return _histograms->getHistogram(name, scheme, specialization);
    }

    EventIdGroupDict<IOHistogram> Hammer::getHistograms(const string& name, const string& scheme,
                                                        const string& specialization) const {
        return _histograms->getHistograms(name, scheme, specialization);
    }

    EventUIDGroup Hammer::getHistogramEventIds(const string& name, const string& scheme,
                                               const string& specialization) const {
        return _histograms->getHistogramEventIds(name, specialization, scheme);
    }

    vector<vector<double>> Hammer::getHistogramBinEdges(const string& name) const {
        return _histograms->getHistogramEdges(name);
    }

    IndexList Hammer::getHistogramShape(const string& name) const {
        return _histograms->getHistogramShape(name);
    }

    bool Hammer::histogramHasUnderOverFlows(const string& name) const {
        return _histograms->getUnderOverFlows(name);
    }

#ifdef HAVE_ROOT

    unique_ptr<TH1D> Hammer::getHistogram1D(const string& name, const string& scheme,
                                            const string& specialization) const {
        return _histograms->getHistogram1D(name, scheme, specialization);
    }

    unique_ptr<TH2D> Hammer::getHistogram2D(const string& name, const string& scheme,
                                            const string& specialization) const {
        return _histograms->getHistogram2D(name, scheme, specialization);
    }

    unique_ptr<TH3D> Hammer::getHistogram3D(const string& name, const string& scheme,
                                            const string& specialization) const {
        return _histograms->getHistogram3D(name, scheme, specialization);
    }

    EventIdGroupDict<unique_ptr<TH1D>> Hammer::getHistograms1D(const string& name, const string& scheme,
                                                               const string& specialization) const {
        return _histograms->getHistograms1D(name, scheme, specialization);
    }

    EventIdGroupDict<unique_ptr<TH2D>> Hammer::getHistograms2D(const string& name, const string& scheme,
                                                               const string& specialization) const {
        return _histograms->getHistograms2D(name, scheme, specialization);
    }

    EventIdGroupDict<unique_ptr<TH3D>> Hammer::getHistograms3D(const string& name, const string& scheme,
                                                               const string& specialization) const {
        return _histograms->getHistograms3D(name, scheme, specialization);
    }

    void Hammer::setHistogram1D(const string& name, const string& scheme, TH1D& histogram,
                                const string& specialization) const {
        _histograms->setHistogram1D(name, scheme, histogram, specialization);
    }

    void Hammer::setHistogram2D(const string& name, const string& scheme, TH2D& histogram,
                                const string& specialization) const {
        _histograms->setHistogram2D(name, scheme, histogram, specialization);
    }

    void Hammer::setHistogram3D(const string& name, const string& scheme, TH3D& histogram,
                                const string& specialization) const {
        _histograms->setHistogram3D(name, scheme, histogram, specialization);
    }

    void Hammer::setHistograms1D(const string& name, const string& scheme,
                                 EventIdGroupDict<unique_ptr<TH1D>>& histograms, const string& specialization) const {
        _histograms->setHistograms1D(name, scheme, histograms, specialization);
    }

    void Hammer::setHistograms2D(const string& name, const string& scheme,
                                 EventIdGroupDict<unique_ptr<TH2D>>& histograms, const string& specialization) const {
        _histograms->setHistograms2D(name, scheme, histograms, specialization);
    }

    void Hammer::setHistograms3D(const string& name, const string& scheme,
                                 EventIdGroupDict<unique_ptr<TH3D>>& histograms, const string& specialization) const {
        _histograms->setHistograms3D(name, scheme, histograms, specialization);
    }

#endif

} // namespace Hammer
