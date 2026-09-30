///
/// @file  ProvidersRepo.cc
/// @brief Interface class for amplitudes, rates, FFs dictionary container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include "Hammer/ProvidersRepo.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/Config/HammerModules.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/Histogram.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {

    // NOLINTBEGIN(readability-function-cognitive-complexity,readability-function-size)
    ProvidersRepo::ProvidersRepo(const SchemeDefinitions* schemeDefs) : _schemeDefs{schemeDefs} {
        ADD_AMPLITUDES(_amplitudes);
        ADD_FORMFACTORS(_formFactors);
        ADD_RATES(_rates);
        PID& pdg = PID::instance();
        _partialWidths = pdg.getPartialWidths();
    }
    // NOLINTEND(readability-function-cognitive-complexity,readability-function-size)

    ProvidersRepo::~ProvidersRepo() noexcept {
        _amplitudes.clear();
        _rates.clear();
        _formFactors.clear();
        _partialWidths.clear();
        _schemeNamesFromPrefixGroup.clear();
        _mWCProviders.clear();
        _mFFErrProviders.clear();
        _schemeDefs = nullptr;
    }

    void ProvidersRepo::defineSettings() {
    }

    void ProvidersRepo::setSettingsHandler(SettingsHandler& sh) {
        SettingsConsumer::setSettingsHandler(sh);
        for (auto& elem : _amplitudes) {
            elem.second->setSettingsHandler(sh);
        }
        for (auto& elem : _rates) {
            elem.second->setSettingsHandler(sh);
        }
        for (auto& elem : _formFactors) {
            for (auto& elem2 : elem.second) {
                elem2->setSettingsHandler(sh);
            }
        }
    }

    void ProvidersRepo::initialize() {
        initializeWCs();
        initializeFFErrs();
    }

    void ProvidersRepo::initializeWCs() {
        for (auto& elem : _amplitudes) {
            auto info = elem.second->getWCInfo();
            if (info.second != IndexLabel::NONE && !info.first.empty()) {
                if (_mWCProviders.find(info.first) == _mWCProviders.end()) {
                    _mWCProviders.insert({info.first, {info.second, elem.second.get()}});
                }
            }
        }
    }

    void ProvidersRepo::initializeFFErrs() {
        for (auto& elem : _formFactors) {
            for (auto& elem2 : elem.second) {
                auto info = elem2->getFFErrInfo();
                if (info.second != IndexLabel::NONE && !info.first.get().empty()) {
                    if (_mFFErrProviders.find(info.first) == _mFFErrProviders.end()) {
                        _mFFErrProviders.insert({info.first, {info.second, elem2.get()}});
                    }
                }
            }
        }
    }

    AmplitudeBase* ProvidersRepo::getAmplitude(PdgId parent, const vector<PdgId>& daughters,
                                               const vector<PdgId>& granddaughters) const {
        auto tmp = combineDaughters(daughters, granddaughters);
        HashId id = processID(parent, tmp);
        auto it = _amplitudes.find(id);
        if (it != _amplitudes.end()) {
            return &(*it->second);
        }
        tmp = combineDaughters(flipSigns(daughters), flipSigns(granddaughters));
        id = processID(flipSign(parent), tmp);
        it = _amplitudes.find(id);
        if (it != _amplitudes.end()) {
            return &(*it->second);
        }
        return nullptr;
    }

    RateBase* ProvidersRepo::getRate(PdgId parent, const vector<PdgId>& daughters,
                                     const vector<PdgId>& granddaughters) const {
        auto tmp = combineDaughters(daughters, granddaughters);
        HashId id = processID(parent, tmp);
        auto it = _rates.find(id);
        if (it != _rates.end()) {
            return &(*it->second);
        }
        tmp = combineDaughters(flipSigns(daughters), flipSigns(granddaughters));
        id = processID(flipSign(parent), tmp);
        it = _rates.find(id);
        if (it != _rates.end()) {
            return &(*it->second);
        }
        return nullptr;
    }

    const double* ProvidersRepo::getPartialWidth(PdgId parent, const vector<PdgId>& daughters,
                                                 const vector<PdgId>& granddaughters) const {
        auto tmp = combineDaughters(daughters, granddaughters);
        HashId id = processID(parent, tmp);
        auto it = _partialWidths.find(id);
        if (it != _partialWidths.end()) {
            return &(it->second);
        }
        tmp = combineDaughters(flipSigns(daughters), flipSigns(granddaughters));
        id = processID(flipSign(parent), tmp);
        it = _partialWidths.find(id);
        if (it != _partialWidths.end()) {
            return &(it->second);
        }
        return nullptr;
    }

    bool ProvidersRepo::isHadronicIdKey(const HashId& processId) const {
        auto ffs = getFormFactor(processId);
        return ffs.size() > 0 ? true : false;
    }
    
    vector<FormFactorBase*> ProvidersRepo::getFormFactor(HashId processId) const {
        vector<FormFactorBase*> res{};
        auto it = _formFactors.find(processId);
        if (it != _formFactors.end()) {
            for (const auto& elem : it->second) {
                res.push_back(&(*elem));
            }
        }
        return res;
    }

    // vector<FormFactorBase*> ProvidersRepo::getFormFactor(PdgId parent, const vector<PdgId>& daughters, const
    // vector<PdgId>& granddaughters) const {
    //     AmplitudeBase* amp = getAmplitude(parent, daughters, granddaughters);
    //     if(amp != nullptr) {
    //         HashId ffId = amp->hadronicId();
    //         return getFormFactor(ffId);
    //     }
    //     else {
    //         return vector<FormFactorBase*>{};
    //     }
    // }

    map<HashId, vector<string>> ProvidersRepo::getFFGroups() const {
        map<HashId, vector<string>> FFAssocs;
        for (const auto& elem : _formFactors) {
            auto res = FFAssocs.insert({elem.first, vector<string>{}});
            if (res.second) {
                for (const auto& elem2 : elem.second) {
                    res.first->second.push_back(elem2->group());
                }
            }
        }
        return FFAssocs;
    }

    map<string, set<string>> ProvidersRepo::getAvailableFFParams() const {
        map<string, set<string>> FFGroupAssocs;
        for (const auto& elem : _formFactors) {
            for (const auto& elem2 : elem.second) {
                auto prefixgroup = elem2->getFFErrPrefixGroup();
                auto itp = FFGroupAssocs.find(prefixgroup.prefix);
                if (itp == FFGroupAssocs.end()) {
                    FFGroupAssocs.insert({prefixgroup.prefix, set<string>{}});
                }
                FFGroupAssocs[prefixgroup.prefix].insert(prefixgroup.group);
            }
        }
        return FFGroupAssocs;
    }

    void ProvidersRepo::showAvailableFFParams() const {
        auto FFGroupAssocs = getAvailableFFParams();
        for (const auto& elem : FFGroupAssocs) {
            string groupset{};
            for (const auto& group : elem.second) {
                groupset += (group + " ");
            }
            MSG_INFO("Known FF parametrizations for '" + elem.first + "': " + groupset);
        }
    }

    void ProvidersRepo::showAvailableFFParams(const string& prefix) const {
        auto FFGroupAssocs = getAvailableFFParams();
        string groupset{};
        auto it = FFGroupAssocs.find(prefix);
        if (it != FFGroupAssocs.end()) {
            for (const auto& group : it->second) {
                groupset += (group + " ");
            }
            MSG_INFO("Known FF parametrizations for '" + prefix + "': " + groupset);
        } else {
            MSG_ERROR("FF class prefix '" + prefix +
                      "' is unknown. Paging Mariner 1 rocket to the courtesy telephone.");
        }
    }

    void ProvidersRepo::initPreSchemeDefs() {
        // this one is required by SchemeDefs init
        auto formFactDups = _schemeDefs->getFFDuplicates();
        SettingsHandler* settings = getSettingsHandler();
        if (settings != nullptr) {
            fixWCValues();
            for (auto& elem : _amplitudes) {
                elem.second->init();
            }
            for (auto& elem : _rates) {
                elem.second->init();
            }
            for (auto& elem : _formFactors) {
                auto itDups = formFactDups.find(elem.first);
                size_t initsize = elem.second.size();
                for (size_t i = 0; i < initsize; ++i) {
                    if (itDups != formFactDups.end()) {
                        auto itList = itDups->second.find(elem.second[i]->group());
                        if (itList != itDups->second.end()) {
                            for (const auto& token : itList->second) {
                                duplicateFormFactor(elem.first, i, token);
                            }
                        }
                    }
                }
                for (auto& elem2 : elem.second) {
                    elem2->init();
                }
            }
            initializeFFErrs(); // so that also the duplicates are initialized
        }
    }

    void ProvidersRepo::initPostSchemeDefs() {
        // this one needs SchemeDefs to be initialized
        _schemeNamesFromPrefixGroup.clear();
        auto oldlevel = Log::getLog("Hammer.SchemeDefinitions").getLevel();
        Log::setLevel("Hammer.SchemeDefinitions", Log::CRITICAL);
        for (auto& elem : _formFactors) {
            auto schemes = _schemeDefs->getFFSchemesForProcess(elem.first);
            for (auto& elem2 : schemes) {
                auto prefgroup = elem.second[elem2.second]->getFFErrPrefixGroup();
                _schemeNamesFromPrefixGroup[prefgroup].insert(elem2.first);
            }
        }
        for (auto& elem : _formFactors) {
            for (auto& elem2 : elem.second) {
                elem2->calcUnits();
            }
        }
        Log::setLevel("Hammer.SchemeDefinitions", oldlevel);
    }

    FFIndex ProvidersRepo::duplicateFormFactor(HashId processId, FFIndex index, const string& newLabel) {
        auto it = _formFactors.find(processId);
        if (it != _formFactors.end()) {
            if (index < it->second.size()) {
                FFIndex newpos = it->second.size();
                it->second.push_back(unique_ptr<FormFactorBase>(it->second[index]->clone(newLabel)));
                return newpos;
            }
        }
        return 0;
    }

    Log& ProvidersRepo::getLog() {
        return Log::getLog("Hammer.ProvidersRepo");
    }

    void ProvidersRepo::fixWCValues() {
        SettingsHandler* settings = getSettingsHandler();
        if (settings != nullptr) {
            for (auto& elem : _amplitudes) {
                string prefix = elem.second->getWCInfo().first;
                auto coeffs = settings->getSettings(prefix);
                auto coeffsNum = settings->getSettings(prefix, WTerm::NUMERATOR);
                auto coeffsDen = settings->getSettings(prefix, WTerm::DENOMINATOR);
                if (coeffs.size() > coeffsNum.size()) {
                    auto it1 = coeffs.begin();
                    auto it2 = coeffsNum.begin();
                    while (it1 != coeffs.end() && it2 != coeffsNum.end()) {
                        if (*it1 < *it2) {
                            settings->cloneSetting(prefix, *it1, *settings->getEntry(prefix, *it1, WTerm::COMMON),
                                                   WTerm::NUMERATOR);
                            settings->resetSetting(prefix, *it1, WTerm::NUMERATOR);
                            ++it1;
                        } else if (*it2 < *it1) {
                            ++it2;
                        } else {
                            ++it1;
                            ++it2;
                        }
                    }
                }
                if (coeffs.size() > coeffsDen.size()) {
                    auto it1 = coeffs.begin();
                    auto it2 = coeffsDen.begin();
                    while (it1 != coeffs.end() && it2 != coeffsDen.end()) {
                        if (*it1 < *it2) {
                            settings->cloneSetting(prefix, *it1, *settings->getEntry(prefix, *it1, WTerm::COMMON),
                                                   WTerm::DENOMINATOR);
                            settings->resetSetting(prefix, *it1, WTerm::DENOMINATOR);
                            ++it1;
                        } else if (*it2 < *it1) {
                            ++it2;
                        } else {
                            ++it1;
                            ++it2;
                        }
                    }
                }
            }
        }
    }

    std::set<std::string> ProvidersRepo::schemeNamesFromPrefixAndGroup(const FFPrefixGroup& value) const {
        return getOrDefault(_schemeNamesFromPrefixGroup, value, {});
    }

    bool ProvidersRepo::checkFFPrefixAndGroup(const FFPrefixGroup& value) const {
        for (const auto& elem : _formFactors) {
            for (const auto& elem2 : elem.second) {
                if (elem2->getFFErrPrefixGroup() == value) {
                    return true;
                }
            }
        }
        return false;
    }

    AmplitudeBase* ProvidersRepo::getWCProvider(const string& wcPrefix) const {
        return getOrDefault(_mWCProviders, wcPrefix, {IndexLabel::NONE, nullptr}).second;
    }

    IndexLabel ProvidersRepo::getWCLabel(const string& wcPrefix) const {
        return getOrDefault(_mWCProviders, wcPrefix, {IndexLabel::NONE, nullptr}).first;
    }

    map<IndexLabel, AmplitudeBase*> ProvidersRepo::getAllWCProviders() const {
        map<IndexLabel, AmplitudeBase*> results;
        for (const auto& elem : _mWCProviders) {
            results.insert(elem.second);
        }
        return results;
    }

    FormFactorBase* ProvidersRepo::getFFErrProvider(const FFPrefixGroup& process) const {
        return getOrDefault(_mFFErrProviders, process, {IndexLabel::NONE, nullptr}).second;
    }

    IndexLabel ProvidersRepo::getFFErrLabel(const FFPrefixGroup& process) const {
        return getOrDefault(_mFFErrProviders, process, {IndexLabel::NONE, nullptr}).first;
    }

    std::map<IndexLabel, SchemeDict<FormFactorBase*>> ProvidersRepo::getAllFFErrProviders() const {
        map<IndexLabel, SchemeDict<FormFactorBase*>> results;
        for (const auto& elem : _mFFErrProviders) {
            auto schemeNames = schemeNamesFromPrefixAndGroup(elem.first);
            for (const auto& elem2 : schemeNames) {
                results[elem.second.first][elem2] = elem.second.second;
            }
        }
        return results;
    }

    bool ProvidersRepo::renameFFEigenvectors(const FFPrefixGroup& process, const vector<string>& names) const {
        string optionPath = process.prefix + process.group;
        SettingsHandler* settings = getSettingsHandler();
        auto* currentValues = getSetting<vector<string>>(optionPath, "ErrNames");
        bool shouldReInitExternalData = false;
        if (currentValues == nullptr) {
            settings->addSetting<vector<string>>(optionPath, "ErrNames", names);
        } else {
            auto* defaultValues =
                settings->getEntry(optionPath, "ErrNames", WTerm::COMMON)->getDefault<vector<string>>();
            if (names.size() > defaultValues->size()) {
                MSG_WARNING("FF Error Names renaming: too many names provided. Excess will be ignored.");
            }
            for (size_t i = 0; i < min(names.size(), defaultValues->size()); ++i) {
                if (!names[i].empty()) {
                    if (names[i] != currentValues->at(i)) {
                        shouldReInitExternalData = true;
                        if (currentValues->at(i) != defaultValues->at(i)) {
                            settings->renameSetting(optionPath, currentValues->at(i), names[i]);
                        } else {
                            settings->addSetting<double>(optionPath, names[i], 0.);
                        }
                        (*currentValues)[i] = names[i];
                    }
                } else {
                    if (currentValues->at(i).empty() && !defaultValues->at(i).empty()) {
                        (*currentValues)[i] = defaultValues->at(i);
                    }
                }
            }
        }
        return shouldReInitExternalData;
    }

} // namespace Hammer
