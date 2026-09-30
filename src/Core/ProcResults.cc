///
/// @file  ProcResults.cc
/// @brief Container for process-related results of weight calculation
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/algorithm/string.hpp>

#include "Hammer/ProcResults.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/SpecializationDefinitions.hh"

using namespace std;

namespace Hammer {

    ProcResults::ProcResults(const Serial::FBProcData* msgreader)
        : _processSquaredAmplitude{}, _processAmplitude{}, _specializations{nullptr} {
        read(msgreader, false);
    }

    // ProcResults::~ProcResults() {
    // }

    const Tensor& ProcResults::processAmplitude(WTerm what) const {
        return _processAmplitude.get(what);
    }

    Tensor& ProcResults::processAmplitude(WTerm what) {
        return _processAmplitude.get(what);
    }

    const Tensor& ProcResults::processAmplitudeSquared(WTerm what) const {
        return _processSquaredAmplitude.get(what);
    }

    Tensor& ProcResults::processAmplitudeSquared(WTerm what) {
        return _processSquaredAmplitude.get(what);
    }

    vector<reference_wrapper<const Tensor>> ProcResults::processFormFactors(const string& name) const {
        vector<reference_wrapper<const Tensor>> result;
        auto it = _processFormFactors.find(name);
        if (it != _processFormFactors.end()) {
            transform(it->second.begin(), it->second.end(), back_inserter(result),
                      [](const Tensor& elem) -> reference_wrapper<const Tensor> { return cref(elem); });
        }
        return result;
    }

    vector<reference_wrapper<Tensor>> ProcResults::processFormFactors(const string& name) {
        vector<reference_wrapper<Tensor>> result;
        auto it = _processFormFactors.find(name);
        if (it != _processFormFactors.end()) {
            transform(it->second.begin(), it->second.end(), back_inserter(result),
                      [](Tensor& elem) -> reference_wrapper<Tensor> { return ref(elem); });
        }
        return result;
    }

    void ProcResults::appendFormFactor(const std::string& schemeName, const Tensor& formfact) {
        _processFormFactors[schemeName].push_back(formfact);
    }

    void ProcResults::appendFormFactor(const std::string& schemeName, Tensor&& formfact) {
        _processFormFactors[schemeName].push_back(std::move(formfact));
    }

    void ProcResults::clearFormFactors() {
        _processFormFactors.clear();
    }

    bool ProcResults::haveFormFactors() const {
        return !_processFormFactors.empty();
    }

    vector<string> ProcResults::availableSchemes() const {
        vector<string> result;
        result.reserve(_processFormFactors.size());
        transform(_processFormFactors.begin(), _processFormFactors.end(), back_inserter(result),
                  [](const pair<string, vector<Tensor>>& val) -> string { return val.first; });
        return result;
    }

    void ProcResults::setSpecDefinitions(SpecializationDefinitions* specs) {
        _specializations = specs;
    }

    const Tensor& ProcResults::processWeight(const string& schemeName, const string& specializationId) const {
        return getOrThrow(getOrThrow(_processWeights, specializationId, RangeError("Invalid specialization Id")),
                          schemeName, RangeError("Invalid scheme name"));
    }

    Tensor& ProcResults::processWeight(const string& schemeName, const string& specializationId) {
        return getOrThrow(getOrThrow(_processWeights, specializationId, RangeError("Invalid specialization Id")),
                          schemeName, RangeError("Invalid scheme name"));
    }

    bool ProcResults::haveWeight(const string& schemeName, const string& specializationId) const {
        return checkExistence(_processWeights, specializationId, schemeName);
    }


    void ProcResults::setProcessWeight(const string& schemeName, const Tensor& weight, const string& specializationId) {
        _processWeights[specializationId][schemeName] = weight;
    }

    void ProcResults::setProcessWeight(const string& schemeName, Tensor&& weight, const string& specializationId) {
        _processWeights[specializationId][schemeName] = std::move(weight);
    }

    void ProcResults::clearWeights() {
        _processWeights.clear();
    }

    void ProcResults::removeSpecializedWeights(const std::string& specializationId) {
        _processWeights.erase(specializationId);
    }

    void ProcResults::removeAllSpecializedWeights() {
        // wish I had C++20 with erase_if...
        for (auto it2 = _processWeights.begin(); it2 != _processWeights.end();) {
            if (it2->first != Spec::none()) {
                it2 = _processWeights.erase(it2);
            } else {
                it2++;
            }
        }
    }

    bool ProcResults::haveWeights() const {
        return !_processWeights.empty();
    }

    void ProcResults::defineSettings() {
        setPath("ProcessWrite");
        addSetting<bool>("Amplitudes", true);
        addSetting<bool>("SquaredAmplitudes", true);
        addSetting<bool>("FormFactors", true);
    }

    void ProcResults::write(flatbuffers::FlatBufferBuilder* msgwriter,
                            flatbuffers::Offset<Serial::FBProcData>* msg) const {
        vector<flatbuffers::Offset<Serial::FBFormFactor>> ffs;
        vector<flatbuffers::Offset<Serial::FBTensor>> wgtsVals;
        vector<string> wgtsNames;
        vector<string> wgtsSpecs;
        ffs.reserve(_processFormFactors.size());
        wgtsVals.reserve(_processWeights.size());
        wgtsNames.reserve(_processWeights.size());
        wgtsSpecs.reserve(_processWeights.size());
        if (isOn("FormFactors")) {
            for (const auto& elem : _processFormFactors) {
                auto name = msgwriter->CreateString(elem.first);
                vector<flatbuffers::Offset<Serial::FBTensor>> ffVec;
                ffVec.reserve(elem.second.size());
                for (const auto& elem2 : elem.second) {
                    flatbuffers::Offset<Serial::FBTensor> val;
                    elem2.write(msgwriter, &val);
                    ffVec.push_back(val);
                }
                auto serialFFVec = msgwriter->CreateVector(ffVec);
                Serial::FBFormFactorBuilder serialFF{*msgwriter};
                serialFF.add_id(name);
                serialFF.add_ffsvals(serialFFVec);
                auto resFF = serialFF.Finish();
                ffs.push_back(resFF);
            }
        }
        for (const auto& elem : _processWeights) {
            if (elem.first == Spec::none() && !isOn("Hammer", "SaveGeneralWeights")) {
                continue; // if CalcGeneralWeights is false, no general weight was added to _processWeights
            }
            if (elem.first != Spec::none() && !_specializations->isUsedInWeights(elem.first)) {
                continue;
            }
            for (const auto& elem2 : elem.second) {
                wgtsNames.push_back(elem2.first);
                wgtsSpecs.push_back(elem.first);
                flatbuffers::Offset<Serial::FBTensor> val;
                elem2.second.write(msgwriter, &val);
                wgtsVals.push_back(val);
            }
        }
        flatbuffers::Offset<Serial::FBTensor> ampnum;
        flatbuffers::Offset<Serial::FBTensor> amp2num;
        flatbuffers::Offset<Serial::FBTensor> ampden;
        flatbuffers::Offset<Serial::FBTensor> amp2den;
        if (isOn("Amplitudes")) {
            _processAmplitude.numerator.write(msgwriter, &ampnum);
            _processAmplitude.denominator.write(msgwriter, &ampden);
        }
        if (isOn("SquaredAmplitudes")) {
            _processSquaredAmplitude.numerator.write(msgwriter, &amp2num);
            _processSquaredAmplitude.denominator.write(msgwriter, &amp2den);
        }
        flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<Serial::FBFormFactor>>> serialFF;
        if (isOn("FormFactors")) {
            serialFF = msgwriter->CreateVector(ffs);
        }
        auto serialWgtN = msgwriter->CreateVectorOfStrings(wgtsNames);
        auto serialWgtV = msgwriter->CreateVector(wgtsVals);
        auto serialWgtS = msgwriter->CreateVectorOfStrings(wgtsSpecs);
        Serial::FBProcDataBuilder serialprocess{*msgwriter};
        if (isOn("Amplitudes")) {
            serialprocess.add_ampnum(ampnum);
            serialprocess.add_ampden(ampden);
        }
        if (isOn("SquaredAmplitudes")) {
            serialprocess.add_amp2num(amp2num);
            serialprocess.add_amp2den(amp2den);
        }
        if (isOn("FormFactors")) {
            serialprocess.add_formfacts(serialFF);
        }
        serialprocess.add_weightames(serialWgtN);
        serialprocess.add_weighvalues(serialWgtV);
        serialprocess.add_weighspecs(serialWgtS);
        *msg = serialprocess.Finish();
    }

    bool ProcResults::read(const Serial::FBProcData* msgreader, bool merge) {
        if (msgreader != nullptr) {
            bool result = true;
            if (!merge && flatbuffers::IsFieldPresent(msgreader, Serial::FBProcData::VT_AMPNUM)) {
                _processAmplitude.numerator.read(msgreader->ampnum());
                _processAmplitude.denominator.read(msgreader->ampden());
            }
            if (!merge && flatbuffers::IsFieldPresent(msgreader, Serial::FBProcData::VT_AMP2NUM)) {
                _processSquaredAmplitude.numerator.read(msgreader->amp2num());
                _processSquaredAmplitude.denominator.read(msgreader->amp2den());
            }
            if (flatbuffers::IsFieldPresent(msgreader, Serial::FBProcData::VT_FORMFACTS)) {
                const auto* formfacts = msgreader->formfacts();
                if (!merge) {
                    _processFormFactors.clear();
                }
                for (unsigned int i = 0; i < formfacts->size() && result; ++i) {
                    const auto* elem = formfacts->Get(i);
                    if (elem->id() == nullptr || elem->ffsvals() == nullptr) {
                        return false;
                    }
                    const char* elemId = elem->id()->c_str();
                    auto it = _processFormFactors.find(elemId);
                    if (it != _processFormFactors.end() && strcmp(elemId, "Denominator") != 0) {
                        MSG_ERROR("Try to merge two process form factor schemes with same name '" + string(elemId) +
                                  "'!");
                        result = false;
                    } else {
                        auto res = _processFormFactors.insert({elemId, vector<Tensor>{}});
                        const auto* ffVec = elem->ffsvals();
                        for (unsigned int j = 0; j < ffVec->size(); ++j) {
                            const auto* elem2 = ffVec->Get(j);
                            Tensor t;
                            t.read(elem2);
                            res.first->second.push_back(std::move(t));
                        }
                    }
                }
            }
            const auto* wName = msgreader->weightames();
            const auto* wVal = msgreader->weighvalues();
            if (wName == nullptr || wVal == nullptr) {
                return false;
            }
            bool hasSpecs = flatbuffers::IsFieldPresent(msgreader, Serial::FBProcData::VT_WEIGHSPECS);
            const auto* wSpcs = hasSpecs ? msgreader->weighspecs() : nullptr;
            if (!merge) {
                _processWeights.clear();
            }
            for (unsigned int i = 0; i < wName->size() && result; ++i) {
                const auto* specStr = wSpcs != nullptr ? wSpcs->Get(i) : nullptr;
                auto id = (specStr != nullptr) ? specStr->c_str() : Spec::none();
                const auto* nameStr = wName->Get(i);
                if (nameStr == nullptr) {
                    return false;
                }
                const auto* name = nameStr->c_str();
                auto it = _processWeights.find(id);
                if (it == _processWeights.end()) {
                    it = _processWeights.insert({id, {}}).first;
                }
                auto it2 = it->second.find(name);
                if (it2 != it->second.end()) {
                    MSG_ERROR("Try to merge two process weights with same name '" + string(name) + "' and same id '" +
                              id + "'!");
                    result = false;
                } else {
                    Tensor t;
                    t.read(wVal->Get(i));
                    it->second.insert({name, t});
                }
            }
            return result;
        }
        return false;
    }

    Log& ProcResults::getLog() {
        return Log::getLog("Hammer.ProcResults");
    }

} // namespace Hammer
