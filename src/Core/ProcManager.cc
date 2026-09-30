///
/// @file  ProcManager.cc
/// @brief Container class for all process related data structures
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/ProcManager.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/RateBase.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/Math/PhaseSpace.hh"
#include "Hammer/ProcRates.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    ProcManager::ProcManager(const Process& inputs)
        : _inputs{new Process{inputs}}, _graph{new ProcGraph{}}, _reqs{new ProcRequirements{}},
          _results{new ProcResults{}}, _processRates{nullptr}, _dictionaries{nullptr} {
        // initSettings();
    }

    ProcManager::ProcManager() : ProcManager{Process{}} {
    }

    ProcManager::ProcManager(const Serial::FBProcess* msgreader) : ProcManager{Process{}} {
        read(msgreader, false);
    }

    ProcManager::ProcManager(ProcManager&& other) noexcept // NOLINT(bugprone-exception-escape)
        : SettingsConsumer{other}, _inputs{other._inputs.release()}, _graph{other._graph.release()},
          _reqs{other._reqs.release()}, _results{other._results.release()}, _processRates{other._processRates},
          _dictionaries{other._dictionaries} {
        // initSettings();
    }

    // ProcManager::~ProcManager() noexcept {
    // }

    void ProcManager::setSettingsHandler(SettingsHandler& sh) {
        SettingsConsumer::setSettingsHandler(sh);
        _results->setSettingsHandler(sh);
    }

    void ProcManager::setSettingsHandler(const SettingsConsumer& sh) {
        SettingsConsumer::setSettingsHandler(sh);
        _results->setSettingsHandler(sh);
    }

    void ProcManager::setDictionary(DictionaryManager* dict) {
        _dictionaries = dict;
    }

    void ProcManager::initialize(DictionaryManager* dictionaries) {
        setDictionary(dictionaries);
        _results->setSpecDefinitions(&_dictionaries->specializationDefs());
        if (_inputs->initialize()) {
            if (_dictionaries->processDefs().isProcessIncluded(_inputs->fullId(), _inputs->getId())) {
                _graph->initialize(dictionaries, _inputs.get());
                auto numAmpls = _reqs->initialize(dictionaries, _inputs.get(), _graph.get());
                if (numAmpls == 0) { //|| _requiredFormFactors.size() == 0 || _requiredRate == nullptr) {
                    _dictionaries->processDefs().forbidProcess(_inputs->getId());
                    _inputs->disable();
                    return;
                }
                for (const auto& id : _reqs->rateIds()) {
                    _processRates = _dictionaries->rates().getProcessRates(id.second);
                    if (_processRates != nullptr && _processRates->empty()) {
                        calcRates(id.first);
                    }
                }
            } else {
                _inputs->disable();
                return;
            }
        }
        /// @todo what else needs to be done here?
        // The 'else'. We'll never get that two hours back.
    }

    void ProcManager::calc() {
        // calc amps & amps^2
        if (isOn("Amplitudes")) {
            pair<bool, bool> first{true, true};
            _results->processAmplitude(WTerm::NUMERATOR) = Tensor{"PurePSNumerator", MD::makeScalar(1.)};
            _results->processAmplitude(WTerm::DENOMINATOR) = Tensor{"PurePSDenominator", MD::makeScalar(1.)};
            for (const auto& elem : _reqs->amplitudes()) {
                auto amplNum = elem.second.amplitudes.numerator;
                auto amplDen = elem.second.amplitudes.denominator;
                auto* numamp = amplNum.ptr;
                auto* denamp = amplDen.ptr;
                if (numamp != nullptr) {
                    auto particles = _graph->getParticleVectors(amplNum.type, elem.first, elem.second.daughterIdx);
                    numamp->eval(get<0>(particles), get<1>(particles), get<2>(particles));
                    if (first.first) {
                        _results->processAmplitude(WTerm::NUMERATOR) = numamp->getTensor();
                        first.first = false;
                    } else {
                        _results->processAmplitude(WTerm::NUMERATOR).dot(numamp->getTensor());
                    }
                    numamp->addRefs();
                }
                if (denamp != nullptr) {
                    if (numamp != denamp) {
                        auto particles = _graph->getParticleVectors(amplDen.type, elem.first, elem.second.daughterIdx);
                        denamp->eval(get<0>(particles), get<1>(particles), get<2>(particles));
                    }
                    if (first.second) {
                        _results->processAmplitude(WTerm::DENOMINATOR) = denamp->getTensor();
                        first.second = false;
                    } else {
                        _results->processAmplitude(WTerm::DENOMINATOR).dot(denamp->getTensor());
                    }
                    denamp->addRefs();
                }
            }
            if (isOn("CheckForNaNs")) {
                if (_results->processAmplitude(WTerm::NUMERATOR).hasNaNs() ||
                    _results->processAmplitude(WTerm::DENOMINATOR).hasNaNs()) {
                    throw NumericalError("NaNs present in an amplitude. Seagulls! (Stop it now)");
                }
            }
        }
        if (isOn("SquaredAmplitudes")) { // && _results->processAmplitude(WTerm::NUMERATOR).rank() > 0) {
            _results->processAmplitudeSquared(WTerm::NUMERATOR) =
                outerSquare(_results->processAmplitude(WTerm::NUMERATOR));
            _results->processAmplitudeSquared(WTerm::NUMERATOR).spinSum();
            _results->processAmplitudeSquared(WTerm::NUMERATOR) *= _reqs->calcCorrectionFactor(WTerm::NUMERATOR);
        }
        if (isOn("SquaredAmplitudes")) { // && _results->processAmplitude(WTerm::DENOMINATOR).rank() > 0) {
            _results->processAmplitudeSquared(WTerm::DENOMINATOR) =
                outerSquare(_results->processAmplitude(WTerm::DENOMINATOR));
            _results->processAmplitudeSquared(WTerm::DENOMINATOR).spinSum();
            _results->processAmplitudeSquared(WTerm::DENOMINATOR) *= _reqs->calcCorrectionFactor(WTerm::DENOMINATOR);
        }
        if (isOn("FormFactors") && _dictionaries != nullptr) {
            // calc ffs
            _results->clearFormFactors();
            const auto& dict = _dictionaries->schemeDefs().getSchemeDefs();
            Tensor noFF{"noFF", MD::makeScalar(1.)};
            for (const auto& elem : _reqs->amplitudes()) {
                auto itFF = _reqs->formFactors().find(elem.first);
                if (itFF != _reqs->formFactors().end()) {
                    AmplitudeBase* ampNum = elem.second.amplitudes.numerator.ptr;
                    AmplitudeBase* ampDen = elem.second.amplitudes.denominator.ptr;
                    HashId id =
                        (ampNum != nullptr)
                            ? ampNum->hadronicId()
                            : ampDen->hadronicId(); // num and den cannot both be nullptr in _requiredFormFactors
                    auto it = _inputs->find(elem.first);
                    if (it == _inputs->end()) {
                        MSG_ERROR("Something really stinky just happened in CalcFormFactors. Can you smell it?");
                        return;
                    }
                    auto daughters = _inputs->getDaughters(it->first, true);
                    const auto& parent = _inputs->getParticle(it->first);
                    auto siblings = _inputs->getSiblings(it->first, true);
                    for (const auto& elemCh : dict) {
                        auto itCh = elemCh.second.find(id);
                        if (itCh != elemCh.second.end()) {
                            auto iFF = itFF->second.find(itCh->second);
                            if (iFF != itFF->second.end()) {
                                if ((elemCh.first == "Denominator" && ampDen == nullptr) ||
                                    (elemCh.first != "Denominator" && ampNum == nullptr)) {
                                    _results->appendFormFactor(elemCh.first, noFF);
                                    continue;
                                }
                                iFF->second->eval(parent, daughters, siblings);
                                _results->appendFormFactor(elemCh.first, iFF->second->getTensor());
                                iFF->second->addRefs();
                            }
                        }
                    }
                } else {
                    for (const auto& elemCh : dict) {
                        _results->appendFormFactor(elemCh.first, noFF);
                    }
                }
            }
        }
        if (isOn("Weights")) {
            // calc denominator
            Tensor denominator = _results->processAmplitudeSquared(WTerm::DENOMINATOR);
            if (_results->haveFormFactors()) {
                if (denominator.hasFFLabels()) {
                    auto ffdens = _results->processFormFactors("Denominator");
                    decltype(ffdens) ffdensvar{};
                    for (auto& elem : ffdens) {
                        if (elem.get().hasFFVarLabels()) { // defer contractions with FF var to last
                            ffdensvar.push_back(elem);
                            continue;
                        }
                        Tensor t = elem;
                        t.outerSquare();
                        denominator.dot(t);
                    }
                    for (const auto& elem : ffdensvar) {
                        Tensor t = elem;
                        t.outerSquare();
                        denominator.dot(t);
                    }
                }
            }
            if (denominator.hasFFVarLabels()) {
                for (const auto& elem : _reqs->denominatorFFEigenVectors()) {
                    denominator.dot(elem);
                }
            }
            if (denominator.hasWCLabels()) {
                for (const auto& elem : _reqs->denominatorWilsonCoeffs()) {
                    denominator.dot(elem);
                }
            }
            if (denominator.rank() > 0) {
                MSG_ERROR("Denominator has uncontracted indices, with rank " + to_string(denominator.rank()) +
                          ". There's a monster in the basement! Check FF input scheme has been correctly set.");
                return;
            }
            double denValue = denominator.element().real();
            // calc weights
            _results->clearWeights();
            auto specNames = _reqs->availableSpecializations();
            // if(specNames.size() == 0){specNames.insert(Spec::none());} //default for general weight calculation
            if (!isOn("Hammer", "CalcGeneralWeights")) {
                specNames.erase(Spec::none());
                if (specNames.empty()) {
                    MSG_ERROR(
                        "CalcGeneralWeights set to false with no WC specializations declared in weights. The hive is "
                        "out of balance.");
                }
            } else {
                specNames.insert(Spec::none()); // default general weights
            }
            for (const auto& specName : specNames) {
                Tensor tbase = _results->processAmplitudeSquared(WTerm::NUMERATOR);
                // Contract with fully specialized WCs, if any.
                vector<reference_wrapper<const Tensor>> specnumsproj{};
                if (tbase.hasWCLabels() && specName != Spec::none()) {
                    auto specnums = _reqs->specializedWilsonCoeffs(specName);
                    for (const auto& elem2 : specnums) {
                        if (elem2.get().rank() == 2) { // defer partial specialization contractions with projectors
                            specnumsproj.push_back(elem2);
                            continue;
                        }
                        tbase.dot(elem2);
                    }
                }
                for (auto& schemeName : _results->availableSchemes()) {
                    auto tsch = tbase;
                    vector<reference_wrapper<Tensor>> ffnumsvar{};
                    // loop for possible multiple FFs in one process
                    if (tsch.hasFFLabels()) {
                        auto ffnums = _results->processFormFactors(schemeName);
                        for (auto& elem2 : ffnums) {
                            if (elem2.get().hasFFVarLabels()) { // defer contractions with FF var to last
                                ffnumsvar.push_back(elem2);
                                continue;
                            }
                            Tensor t = elem2;
                            t.outerSquare();
                            tsch.dot(t);
                        }
                    }
                    for (const auto& elem2 : specnumsproj) {
                        tsch.dot(elem2);
                    }
                    for (auto& elem2 : ffnumsvar) {
                        Tensor t = elem2;
                        t.outerSquare();
                        tsch.dot(t);
                    }
                    tsch *= 1. / denValue;
                    _results->setProcessWeight(schemeName, tsch, specName);
                }
            }
        }
    }

    const ProcResults& ProcManager::results() const {
        return *_results;
    }

    ProcResults& ProcManager::modifiableRsults() {
        return *_results;
    }

    const Process& ProcManager::inputs() const {
        return *_inputs;
    }

    void ProcManager::defineSettings() {
        setPath("ProcessCalc");
        addSetting<bool>("FormFactors", true);
        addSetting<bool>("SquaredAmplitudes", true);
        addSetting<bool>("Amplitudes", true);
        addSetting<bool>("Rates", true);
        addSetting<bool>("Weights", true);
        addSetting<bool>("CheckForNaNs", false);
    }

    void ProcManager::write(flatbuffers::FlatBufferBuilder* msgwriter,
                            flatbuffers::Offset<Serial::FBProcess>* msg) const {
        flatbuffers::Offset<Serial::FBProcIDs> serialids;
        _inputs->write(msgwriter, &serialids);
        flatbuffers::Offset<Serial::FBProcData> serialdata;
        _results->write(msgwriter, &serialdata);
        Serial::FBProcessBuilder serialprocess{*msgwriter};
        serialprocess.add_ids(serialids);
        serialprocess.add_data(serialdata);
        *msg = serialprocess.Finish();
    }

    bool ProcManager::read(const Serial::FBProcess* msgreader, bool merge) {
        if (msgreader != nullptr) {
            bool res = true;
            if (!merge) {
                _inputs->read(msgreader->ids());
            }
            if (_dictionaries == nullptr) {
                return false;
            }
            if (_dictionaries->processDefs().isProcessIncluded(_inputs->fullId(), _inputs->getId())) {
                res &= _results->read(msgreader->data(), merge);
            }
            return res;
        }
        return false;
    }

    Tensor ProcManager::calcPSTensor(const ParticleIndex& parent) const {
        if (!_inputs->isParent(parent)) {
            MSG_ERROR("Something really stinky just happened in calcPSTensor. Can you smell it?");
        }
        Tensor res{"PS", MD::makeEmptyScalar()};
        const PID& pdg = PID::instance();
        auto pmassPDG = pdg.getMass(_inputs->getParticle(parent).pdgId());
        auto daughters = _inputs->getDaughters(parent);
        vector<double> dmassesPDG;
        for (auto& elem : daughters) {
            dmassesPDG.push_back(pdg.getMass(elem.pdgId()));
        }
        res.element({}) =
            0.5 * phaseSpaceNBody(pmassPDG, dmassesPDG) * pow(pmassPDG, 5 - 2 * static_cast<double>(daughters.size()));
        return res;
    }

    Tensor ProcManager::calcRateTensor(const ParticleIndex& parent, const map<HashId, FFIndex>& ffmaps,
                                       const string& scheme) const {
        auto itr = _reqs->rates().find(parent);
        if (itr != _reqs->rates().end()) {
            RateBase* rateptr = itr->second;
            rateptr->calcTensor();
            Tensor rateT = rateptr->getTensor();
            auto itf = _reqs->formFactors().find(parent);
            auto itidx = ffmaps.find(rateptr->hadronicId());
            if (itf != _reqs->formFactors().end() && itidx != ffmaps.end()) {
                auto itFF = itf->second.find(itidx->second);
                const auto& points = rateptr->getEvaluationPoints();
                Tensor rateFF = itFF->second->getFFPSIntegrand(points);
                rateT.dot(rateFF);
            }
            if (rateT.hasFFLabels()) {
                MSG_ERROR("Unable to contract rate tensor for vertex with parent " +
                          to_string(_inputs->getParticle(parent).pdgId()) + " with scheme " + scheme +
                          ". Eject! Eject!");
                return Tensor{"EmptyRateTensor", MD::makeEmptyScalar()};
            }
            return rateT;
        }
        auto itpw = _reqs->partialWidths().find(parent);
        if (itpw != _reqs->partialWidths().end()) {
            Tensor pw{"PW", MD::makeEmptyScalar()};
            pw.element({}) = *(itpw->second);
            return pw;
        }
        MSG_ERROR("Unable to compute rate tensor for vertex with parent " +
                  to_string(_inputs->getParticle(parent).pdgId()) + " and scheme " + scheme + ". Pull up! Pull up!");
        return Tensor{"EmptyRateTensor", MD::makeEmptyScalar()};
    }

    Tensor ProcManager::calcSpecRateTensor(const string& spec, const string& scheme, const Tensor& rateten) const {
        if (rateten.hasWCLabels() && spec != Spec::none() && scheme != "Denominator") {
            auto res = rateten;
            auto specnums = _reqs->specializedWilsonCoeffs(spec);
            for (const auto& elem : specnums) {
                res.dot(elem);
            }
            return res;
        }
        return rateten;
    }

    void ProcManager::calcRates(const ParticleIndex& parent) const {
        if (!isOn("Rates")) {
            return;
        }
        auto it = _reqs->amplitudes().find(parent);
        if (it == _reqs->amplitudes().end()) {
            return;
        } // Daughter vertices to be filled with parent
        auto amplNum = (it->second).amplitudes.numerator;
        auto amplDen = (it->second).amplitudes.denominator;
        auto daughter = (it->second).daughterIdx;
        // Specializations (same logic as in calc())
        auto specsNum = _reqs->availableSpecializations();
        if (!isOn("Hammer", "GeneralRates")) {
            specsNum.erase(Spec::none());
            if (specsNum.empty()) {
                MSG_ERROR(
                    "GeneralRates set to false with no WC specializations declared in weights. The hive is out of "
                    "balance.");
            }
        } else {
            specsNum.insert(Spec::none()); // default general rate
        }
        set<string> specsDen{Spec::none()}; // for denominator general rate only
        // FF Schemes (kept as outer loop: only do each FF integration once!)
        SchemeDict<map<HashId, FFIndex>> schemes = _dictionaries->schemeDefs().getSchemeDefs();
        for (auto& scheme : schemes) {
            auto ampl = (scheme.first == "Denominator") ? amplDen : amplNum;
            auto specs = (scheme.first == "Denominator") ? specsDen : specsNum;
            switch (ampl.type) {
            case AmplType::VERTEX: {
                if (ampl.ptr == nullptr) {
                    auto rateten = calcPSTensor(parent);
                    for (const auto& spec : specs) {
                        (*_processRates)[spec][scheme.first] = rateten;
                    }
                    //_processRates->insert({scheme.first, calcPSTensor(parent)});
                } else {
                    auto rateten = calcRateTensor(parent, scheme.second, scheme.first);
                    for (const auto& spec : specs) {
                        (*_processRates)[spec][scheme.first] = calcSpecRateTensor(spec, scheme.first, rateten);
                    }
                    //_processRates->insert({scheme.first, calcRateTensor(parent, scheme.second, scheme.first)});
                }
                break;
            }
            case AmplType::FULLEDGE: {
                auto* dprocessRates = _dictionaries->rates().getProcessRates(_reqs->rateIds().at(daughter));
                if (ampl.ptr == nullptr) {
                    auto rateten = calcPSTensor(parent);
                    for (const auto& spec : specs) {
                        (*_processRates)[spec][scheme.first] = rateten;
                    }
                    auto drateten = calcPSTensor(daughter);
                    for (const auto& spec : specs) {
                        (*dprocessRates)[spec][scheme.first] = drateten;
                    }
                    //_processRates->insert({scheme.first, calcPSTensor(parent)});
                    // dprocessRates->insert({scheme.first, calcPSTensor(daughter)});
                } else {
                    auto rateten = calcRateTensor(parent, scheme.second, scheme.first);
                    for (const auto& spec : specs) {
                        (*_processRates)[spec][scheme.first] = calcSpecRateTensor(spec, scheme.first, rateten);
                    }
                    auto drateten = calcRateTensor(daughter, scheme.second, scheme.first);
                    for (const auto& spec : specs) {
                        (*dprocessRates)[spec][scheme.first] = calcSpecRateTensor(spec, scheme.first, drateten);
                    }
                    //_processRates->insert({scheme.first, calcRateTensor(parent, scheme.second, scheme.first)});
                    // dprocessRates->insert({scheme.first, calcRateTensor(daughter, scheme.second, scheme.first)});
                }
                break;
            }
            case AmplType::PARENTEDGE: {
                auto* dprocessRates = _dictionaries->rates().getProcessRates(_reqs->rateIds().at(daughter));
                auto rateten = calcRateTensor(parent, scheme.second, scheme.first);
                for (const auto& spec : specs) {
                    (*_processRates)[spec][scheme.first] = calcSpecRateTensor(spec, scheme.first, rateten);
                }
                auto drateten = calcPSTensor(daughter);
                for (const auto& spec : specs) {
                    (*dprocessRates)[spec][scheme.first] = drateten;
                }
                //_processRates->insert({scheme.first, calcRateTensor(parent, scheme.second, scheme.first)});
                // dprocessRates->insert({scheme.first, calcPSTensor(daughter)});
                break;
            }
            case AmplType::DAUGHTEREDGE: {
                auto* dprocessRates = _dictionaries->rates().getProcessRates(_reqs->rateIds().at(daughter));
                auto rateten = calcPSTensor(parent);
                for (const auto& spec : specs) {
                    (*_processRates)[spec][scheme.first] = rateten;
                }
                auto drateten = calcRateTensor(daughter, scheme.second, scheme.first);
                for (const auto& spec : specs) {
                    (*dprocessRates)[spec][scheme.first] = calcSpecRateTensor(spec, scheme.first, drateten);
                }
                //_processRates->insert({scheme.first, calcPSTensor(parent)});
                // dprocessRates->insert({scheme.first, calcRateTensor(daughter, scheme.second, scheme.first)});
            }
            }
        }
    }

    //    NumDenPair<double> ProcessManager::getPSRates() const {
    //        return _processPhaseSpaceRates;
    //    }

    Log& ProcManager::getLog() {
        return Log::getLog("Hammer.ProcManager");
    }

} // namespace Hammer
