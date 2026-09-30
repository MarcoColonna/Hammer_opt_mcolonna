///
/// @file  ExternalData.cc
/// @brief Container class for values of WC and FF vectors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <algorithm>

#include "Hammer/ExternalData.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Math/Histogram.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/SpecializationDefinitions.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    ExternalData::ExternalData(const IWCFFErrProviders* provs, const SpecializationDefinitions* specs)
        : _providers{provs}, _specializations{specs} {
        _processWilsonCoefficients.reset(new processWilsonCoefficientsType{});
        _processFFEigenVectors.reset(new processFFEigenVectorsType{});
        _externalVectors.reset(new externalVectorsType{});
    }

    ExternalData::~ExternalData() noexcept {
        if (_processWilsonCoefficients.get() != nullptr) {
            _processWilsonCoefficients.reset();
        }
        if (_processFFEigenVectors.get() != nullptr) {
            _processFFEigenVectors.reset();
        }
        if (_externalVectors.get() != nullptr) {
            _externalVectors.reset();
        }
        _schemes.clear();
    }

    void ExternalData::defineSettings() {
    }

    const Tensor& ExternalData::getExternalVectors(const string& schemeName, LabelsList labels,
                                                   const string& specializationId) const {
        // sort(labels.begin(), labels.end(), greater<IndexLabel>());
        if (_externalVectors.get() == nullptr) {
            _externalVectors.reset(new externalVectorsType{});
        }
        auto& ext = (*_externalVectors)[specializationId][schemeName];
        auto ite = ext.find(labels);
        if (ite != ext.end()) {
            return ite->second;
        }
        if (labels.empty()) {
            auto res = ext.insert({labels, Tensor{"External", MD::makeScalar(1.)}});
            return res.first->second;
        }
        vector<pair<MD::SharedTensorData, bool>> elements;
        elements.reserve(labels.size());
        size_t index = (schemeName == "Denominator") ? 1 : 0;
        for (auto elem : labels) {
            bool hc = elem < 0;
            IndexLabel lab = (elem > 0) ? elem : static_cast<IndexLabel>(-elem);
            if (_processWilsonCoefficients.get() != nullptr) {
                auto it = _processWilsonCoefficients->find(lab);
                if (it != _processWilsonCoefficients->end()) {
                    elements.emplace_back(it->second[index], hc);
                    continue;
                }
            }
            if (_processFFEigenVectors.get() != nullptr) {
                auto it2 = _processFFEigenVectors->find(lab);
                if (it2 != _processFFEigenVectors->end()) {
                    auto it3 = it2->second.find(schemeName);
                    if (it3 != it2->second.end()) {
                        elements.emplace_back(it3->second, hc);
                        continue;
                    }
                }
            }
        }
        auto res = ext.emplace(labels, Tensor{"External", std::move(elements)});
        return res.first->second;
    }

    MD::SharedTensorData ExternalData::getWilsonCoefficients(PdgId parent, const vector<PdgId>& daughters,
                                                             const vector<PdgId>& granddaughters, WTerm what) const {
        const AmplitudeBase* amp = _providers->getAmplitude(parent, daughters, granddaughters);
        return getWilsonCoefficients(amp, what);
    }

    MD::SharedTensorData ExternalData::getWilsonCoefficients(const AmplitudeBase* ampl, WTerm what) const {
        if (_processWilsonCoefficients.get() != nullptr) {
            auto it = _processWilsonCoefficients->find(ampl->getWCInfo().second);
            if (it != _processWilsonCoefficients->end()) {
                switch (what) {
                case WTerm::NUMERATOR:
                    return it->second[0];
                case WTerm::DENOMINATOR:
                    return it->second[1];
                case WTerm::COMMON:
                    throw Error("Invalid option");
                }
            }
        }
        return nullptr;
    }

    void ExternalData::setWilsonCoefficients(const string& prefixName, const vector<complex<double>>& values,
                                             WTerm what, bool updateSettings) {
        if (!updateSettings && what != WTerm::NUMERATOR) {
            MSG_ERROR("Only numerator Wilson Coefficients can be updated locally! Skipping");
            return;
        }
        if (_processWilsonCoefficients.get() == nullptr) {
            initWilsonCoefficients();
        }
        auto* ampl = _providers->getWCProvider(prefixName);
        if (ampl != nullptr) {
            auto label = _providers->getWCLabel(prefixName);
            auto itWC = _processWilsonCoefficients->find(label);
            ASSERT(itWC != _processWilsonCoefficients->end());
            if (updateSettings) {
                ampl->updateWCSettings(values, what);
            }
            switch (what) {
            case WTerm::NUMERATOR:
                ampl->updateWCTensor(values, itWC->second[0]);
                break;
            case WTerm::DENOMINATOR:
                ampl->updateWCTensor(values, itWC->second[1]);
                break;
            case WTerm::COMMON:
                ampl->updateWCTensor(values, itWC->second[0]);
                ampl->updateWCTensor(values, itWC->second[1]);
                break;
            }
        } else {
            auto prefixId = SpecPrefixId::toSpecPrefixId(prefixName);
            if (prefixId.hasSpecId()) {
                ASSERT(what == WTerm::NUMERATOR); // specialization only valid in numerator
                auto& wcsub = const_cast<WCSpecialization&>(_specializations->getSpecialization(prefixId));
                auto label = wcsub.getFullLabel();
                auto itWC = _processWilsonCoefficients->find(label);
                ASSERT(itWC != _processWilsonCoefficients->end());
                vector<complex<double>> vec;
                if (updateSettings) {
                    wcsub.updateWCSpecSettings(values);
                    vec = wcsub.getWCSpecVectorFromSettings();
                } else {
                    vec.reserve(values.size() + 1);
                    vec.emplace_back(1.);
                    copy(values.begin(), values.end(), back_inserter(vec));
                }
                wcsub.updateWCSpecTensor(vec, itWC->second[0]);
            } else {
                // TODO: fix this
            }
        }
    }

    void ExternalData::setWilsonCoefficients(const string& prefixName, const map<string, complex<double>>& values,
                                             WTerm what, bool updateSettings) {
        if (!updateSettings && what != WTerm::NUMERATOR) {
            MSG_ERROR("Only numerator Wilson Coefficients can be updated locally! Skipping");
            return;
        }
        if (_processWilsonCoefficients.get() == nullptr) {
            initWilsonCoefficients();
        }
        auto* ampl = _providers->getWCProvider(prefixName);
        if (ampl != nullptr) {
            auto label = _providers->getWCLabel(prefixName);
            auto itWC = _processWilsonCoefficients->find(label);
            ASSERT(itWC != _processWilsonCoefficients->end());
            vector<complex<double>> vec;
            if (updateSettings) {
                ampl->updateWCSettings(values, what);
                vec = ampl->getWCVectorFromSettings(what);
            } else {
                if (values.find("SM") == values.end()) {
                    auto tempValues = values;
                    tempValues["SM"] = 1.0;
                    vec = ampl->getWCVectorFromDict(tempValues);
                } else {
                    vec = ampl->getWCVectorFromDict(values);
                }
            }
            switch (what) {
            case WTerm::NUMERATOR:
                ampl->updateWCTensor(vec, itWC->second[0]);
                break;
            case WTerm::DENOMINATOR:
                ampl->updateWCTensor(vec, itWC->second[1]);
                break;
            case WTerm::COMMON:
                ampl->updateWCTensor(vec, itWC->second[0]);
                ampl->updateWCTensor(vec, itWC->second[1]);
                break;
            }
        } else {
            auto prefixId = SpecPrefixId::toSpecPrefixId(prefixName);
            if (prefixId.hasSpecId()) {
                ASSERT(what == WTerm::NUMERATOR); // specialization only valid in numerator
                auto& wcsub = const_cast<WCSpecialization&>(_specializations->getSpecialization(prefixId));
                auto label = wcsub.getFullLabel();
                auto itWC = _processWilsonCoefficients->find(label);
                ASSERT(itWC != _processWilsonCoefficients->end());
                vector<complex<double>> vec;
                if (updateSettings) {
                    wcsub.updateWCSpecSettings(values);
                    vec = wcsub.getWCSpecVectorFromSettings();
                } else {
                    vec = wcsub.getWCSpecVectorFromDict(values);
                }
                wcsub.updateWCSpecTensor(vec, itWC->second[0]);
            } else {
                // TODO: fix this
            }
        }
    }

    void ExternalData::resetWilsonCoefficients(const string& prefixName, WTerm what) {
        auto* ampl = _providers->getWCProvider(prefixName);
        if (ampl != nullptr) {
            auto label = _providers->getWCLabel(prefixName);
            ASSERT(_processWilsonCoefficients.get() != nullptr);
            auto itWC = _processWilsonCoefficients->find(label);
            ASSERT(itWC != _processWilsonCoefficients->end());
            auto vec = ampl->resetWCSettings(what);
            switch (what) {
            case WTerm::NUMERATOR:
                ampl->updateWCTensor(vec, itWC->second[0]);
                break;
            case WTerm::DENOMINATOR:
                ampl->updateWCTensor(vec, itWC->second[1]);
                break;
            case WTerm::COMMON:
                ampl->updateWCTensor(vec, itWC->second[0]);
                ampl->updateWCTensor(vec, itWC->second[1]);
                break;
            }
        } else {
            auto prefixId = SpecPrefixId::toSpecPrefixId(prefixName);
            if (prefixId.hasSpecId()) {
                ASSERT(what == WTerm::NUMERATOR); // sub-specialization only valid in numerator
                auto& wcsub = const_cast<WCSpecialization&>(_specializations->getSpecialization(prefixId));
                auto label = wcsub.getFullLabel();
                ASSERT(_processWilsonCoefficients.get() != nullptr);
                auto itWC = _processWilsonCoefficients->find(label);
                ASSERT(itWC != _processWilsonCoefficients->end());
                auto vec = wcsub.resetWCSpecialization();
                wcsub.updateWCSpecTensor(vec, itWC->second[0]);
            } else {
            }
        }
    }

    map<string, complex<double>> ExternalData::retrieveWilsonCoefficients(const string& prefixName, WTerm what) const {
        auto* ampl = _providers->getWCProvider(prefixName);
        if (ampl != nullptr) {
            return ampl->retrieveWCSettings(what);
        }
        auto prefixId = SpecPrefixId::toSpecPrefixId(prefixName);
        if (prefixId.hasSpecId()) {
            ASSERT(what == WTerm::NUMERATOR); // sub-specialization only valid in numerator
            auto& wcsub = const_cast<WCSpecialization&>(_specializations->getSpecialization(prefixId));
            return wcsub.retrieveWCSpecialization();
        }
        MSG_WARNING(""); /// TODO: finish this
        return {};
    }

    MD::SharedTensorData ExternalData::getFFEigenVectors(FormFactorBase* ff, const string& schemeName) const {
        ASSERT(_processFFEigenVectors.get() != nullptr);
        auto it = _processFFEigenVectors->find(ff->getFFErrInfo().second);
        if (it != _processFFEigenVectors->end()) {
            auto it2 = it->second.find(schemeName);
            if (it2 != it->second.end()) {
                return it2->second;
            }
        }
        return nullptr;
    }


    void ExternalData::setFFEigenVectors(const FFPrefixGroup& process, const vector<double>& values,
                                         bool updateSettings) {
        auto* ff = _providers->getFFErrProvider(process);
        if (ff != nullptr) {
            if (_processFFEigenVectors.get() == nullptr) {
                initFormFactorErrors();
            }
            auto label = _providers->getFFErrLabel(process);
            auto itFFE1 = _processFFEigenVectors->find(label);
            ASSERT(itFFE1 != _processFFEigenVectors->end());
            if (updateSettings) {
                ff->updateFFErrSettings(values);
            }
            auto schemes = _providers->schemeNamesFromPrefixAndGroup(process);
            vector<double> tmpvalues;
            tmpvalues.reserve(values.size() + 1);
            tmpvalues.push_back(1.);
            tmpvalues.insert(tmpvalues.end(), values.begin(), values.end());
            for (const auto& elem : schemes) {
                auto itFFE2 = itFFE1->second.find(elem);
                if (itFFE2 != itFFE1->second.end()) {
                    ff->updateFFErrTensor(tmpvalues, itFFE2->second);
                }
            }
        }
    }

    void ExternalData::setFFEigenVectors(const FFPrefixGroup& process, const map<string, double>& values,
                                         bool updateSettings) {
        auto* ff = _providers->getFFErrProvider(process);
        if (ff != nullptr) {
            if (_processFFEigenVectors.get() == nullptr) {
                initFormFactorErrors();
            }
            auto label = _providers->getFFErrLabel(process);
            auto itFFE1 = _processFFEigenVectors->find(label);
            ASSERT(itFFE1 != _processFFEigenVectors->end());
            vector<double> vec;
            if (updateSettings) {
                ff->updateFFErrSettings(values);
                vec = ff->getErrVectorFromSettings();
            } else {
                vec = ff->getErrVectorFromDict(values);
            }
            auto schemes = _providers->schemeNamesFromPrefixAndGroup(process);
            for (const auto& elem : schemes) {
                auto itFFE2 = itFFE1->second.find(elem);
                if (itFFE2 != itFFE1->second.end()) {
                    ff->updateFFErrTensor(vec, itFFE2->second);
                }
            }
        }
    }


    void ExternalData::resetFFEigenVectors(const FFPrefixGroup& process) {
        auto* ff = _providers->getFFErrProvider(process);
        if (ff != nullptr) {
            auto label = _providers->getFFErrLabel(process);
            ASSERT(_processFFEigenVectors.get() != nullptr);
            auto itFFE1 = _processFFEigenVectors->find(label);
            ASSERT(itFFE1 != _processFFEigenVectors->end());
            auto n = ff->getErrVectorFromSettings().size();
            vector<double> tmpvalues(n - 1);
            tmpvalues.reserve(n);
            ff->updateFFErrSettings(tmpvalues);
            tmpvalues.insert(tmpvalues.begin(), 1.);
            auto schemes = _providers->schemeNamesFromPrefixAndGroup(process);
            for (const auto& elem : schemes) {
                auto itFFE2 = itFFE1->second.find(elem);
                if (itFFE2 != itFFE1->second.end()) {
                    ff->updateFFErrTensor(tmpvalues, itFFE2->second);
                }
            }
        }
    }

    map<string, double> ExternalData::retrieveFFEigenvectors(const FFPrefixGroup& process) const {
        auto* ff = _providers->getFFErrProvider(process);
        if (ff != nullptr) {
            return ff->retrieveFFErrSettings();
        }
        MSG_ERROR("Process class prefix and group '" + process.get() +
                  "' not found in retrieving FFEigenvectors values. Stormbreaker!");
        return {};
    }

    MD::SharedTensorData ExternalData::getTempFFEigenVectors(const FFPrefixGroup& process,
                                                             const vector<double>& values) const {
        auto* ff = _providers->getFFErrProvider(process);
        MD::SharedTensorData t{new MD::ScalarContainer{}};
        if (ff != nullptr) {
            ff->updateFFErrTensor(values, t);
        }
        return t;
    }

    MD::SharedTensorData ExternalData::getTempFFEigenVectors(const FFPrefixGroup& process,
                                                             const map<string, double>& values) const {
        auto* ff = _providers->getFFErrProvider(process);
        MD::SharedTensorData t{new MD::ScalarContainer{}};
        if (ff != nullptr) {
            auto vec = ff->getErrVectorFromDict(values);
            ff->updateFFErrTensor(vec, t);
        }
        return t;
    }


    void ExternalData::initWilsonCoefficients() {
        _processWilsonCoefficients.reset(new processWilsonCoefficientsType{});
        // _processSpecializedWilsonCoefficients->clear();
        for (auto& elem : _providers->getAllWCProviders()) {
            auto res = _processWilsonCoefficients->emplace(elem.first, array<MD::SharedTensorData, 2>{});
            if (res.second) {
                auto vecN = elem.second->getWCVectorFromSettings(WTerm::NUMERATOR);
                elem.second->updateWCTensor(vecN, res.first->second[0]);
                auto vecD = elem.second->getWCVectorFromSettings(WTerm::DENOMINATOR);
                elem.second->updateWCTensor(vecD, res.first->second[1]);
            }
        }
        for (const auto& elem : _specializations->specializationIds()) {
            for (const auto& elem2 : _specializations->getSpecializations(elem)) {
                if (elem2.second.isPartialSpecialization()) {
                    auto res = _processWilsonCoefficients->emplace(elem2.second.getFullLabel(),
                                                                   array<MD::SharedTensorData, 2>{});
                    if (res.second) {
                        auto vecN = elem2.second.getWCSpecVectorFromSettings();
                        elem2.second.updateWCSpecTensor(vecN, res.first->second[0]);
                    }
                }
            }
        }
    }

    void ExternalData::initFormFactorErrors() {
        _processFFEigenVectors.reset(new processFFEigenVectorsType{});
        for (auto& elem : _providers->getAllFFErrProviders()) {
            auto res = _processFFEigenVectors->emplace(elem.first, SchemeDict<MD::SharedTensorData>{});
            if (res.second) {
                for (auto& elem2 : elem.second) {
                    auto res2 = res.first->second.emplace(elem2.first, MD::SharedTensorData{});
                    if (res2.second) {
                        bool useDefault = (elem2.first == "Denominator");
                        elem2.second->bindErrNames(); // this is to make sure errnames are valid during reload.
                        auto vecN = elem2.second->getErrVectorFromSettings(useDefault);
                        elem2.second->updateFFErrTensor(vecN, res2.first->second);
                    }
                }
            }
        }
    }

    void ExternalData::reInitFormFactorErrors() {
        initFormFactorErrors();
    }

    void ExternalData::initExternalVectors() {
        bool shouldCalcGeneral = isOn("Hammer", "CalcGeneralWeights");
        _externalVectors.reset(new externalVectorsType{});
        auto ids = _specializations->specializationIds();
        if (ids.empty()) {
            ids.insert(Spec::none());
        } else if (!shouldCalcGeneral) {
            ids.erase(Spec::none());
        }
        _externalVectors->reserve(ids.size());
        for (const auto& elem : ids) {
            auto it = _externalVectors->insert({elem, SpecializationDict<UMap<LabelsList, Tensor>>{}});
            it.first->second.reserve(_schemes.size() + 1);
            for (const auto& elem2 : _schemes) {
                it.first->second.insert({elem2, UMap<LabelsList, Tensor>{}});
            }
            it.first->second.insert({"Denominator", UMap<LabelsList, Tensor>{}});
        }
    }

    void ExternalData::init(vector<string> schemeNames) {
        SettingsHandler* settings = getSettingsHandler();
        _schemes = std::move(schemeNames);
        initExternalVectors();
        if (settings != nullptr) {
            initWilsonCoefficients();
            initFormFactorErrors();
        }
    }

    Log& ExternalData::getLog() {
        return Log::getLog("Hammer.ExternalData");
    }


} // namespace Hammer
