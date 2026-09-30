///
/// @file  SchemeDefinitions.cc
/// @brief Container class for Scheme Definitions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/algorithm/string.hpp>

#include "Hammer/Tools/Logging.hh"
#include "yaml-cpp/yaml.h"

#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {


    Log& SchemeDefinitions::getLog() {
        return Log::getLog("Hammer.SchemeDefinitions");
    }

    void SchemeDefinitions::write(flatbuffers::FlatBufferBuilder* msgwriter,
                                  vector<flatbuffers::Offset<Serial::FBFFScheme>>* schemes) const {
        schemes->reserve(_formFactorSchemes.size());
        for (const auto& elem : _formFactorSchemes) {
            vector<uint16_t> chs;
            vector<uint64_t> ids;
            vector<string> keys;
            vector<string> vals;
            for (const auto& elem2 : elem.second) {
                chs.push_back(static_cast<uint16_t>(elem2.second));
                ids.push_back(elem2.first);
            }
            if (elem.first == "Denominator") {
                for (const auto& elem2 : _formFactorBase) {
                    keys.push_back(elem2.first);
                    vals.push_back(elem2.second);
                }
            } else {
                auto it = _formFactorSchemeNames.find(elem.first);
                for (const auto& elem2 : it->second) {
                    keys.push_back(elem2.first);
                    vals.push_back(elem2.second);
                }
            }
            auto serialChs = msgwriter->CreateVector(chs);
            auto serialIds = msgwriter->CreateVector(ids);
            auto serialProcs = msgwriter->CreateVectorOfStrings(keys);
            auto serialGroups = msgwriter->CreateVectorOfStrings(vals);
            auto name = msgwriter->CreateString(elem.first);
            Serial::FBFFSchemeBuilder serialFF{*msgwriter};
            serialFF.add_name(name);
            serialFF.add_ffids(serialIds);
            serialFF.add_ffschemes(serialChs);
            serialFF.add_ffproc(serialProcs);
            serialFF.add_ffgroup(serialGroups);
            auto resFF = serialFF.Finish();
            schemes->push_back(resFF);
        }
    }

    bool SchemeDefinitions::read(const Serial::FBHeader* msgreader, bool merge) {
        const auto* schemes = msgreader->ffschemes();
        if (schemes == nullptr) {
            return false;
        }
        if (!merge) {
            _formFactorSchemes.clear();
            _formFactorSchemeNames.clear();
            _formFactorBase.clear();
        }
        for (unsigned int i = 0; i < schemes->size(); ++i) {
            if (schemes->Get(i)->name() == nullptr) {
                return false;
            }
            string name = schemes->Get(i)->name()->c_str();
            const auto* procs = schemes->Get(i)->ffproc();
            const auto* groups = schemes->Get(i)->ffgroup();
            if (procs == nullptr || groups == nullptr) {
                return false;
            }
            bool result = true;
            if (name == "Denominator") {
                uint8_t denDifferentCount = 0;
                for (unsigned int j = 0; j < procs->size(); ++j) {
                    const auto* pStr = procs->Get(j);
                    const auto* gStr = groups->Get(j);
                    if (pStr == nullptr || gStr == nullptr) {
                        return false;
                    }
                    auto itd = _formFactorBase.find(pStr->c_str());
                    if (itd != _formFactorBase.end()) {
                        if (itd->second != gStr->c_str()) {
                            if (denDifferentCount < 5) {
                                MSG_WARNING("Scheme definition for 'Denominator' differs.");
                            }
                            ++denDifferentCount;
                        }
                    }
                    _formFactorBase.emplace(pStr->c_str(), gStr->c_str());
                }
            } else {
                auto itsn = _formFactorSchemeNames.find(name);
                if (itsn == _formFactorSchemeNames.end()) {
                    itsn = _formFactorSchemeNames.emplace(name, map<string, string>{}).first;
                }
                for (unsigned int j = 0; j < procs->size(); ++j) {
                    const auto* pStr = procs->Get(j);
                    const auto* gStr = groups->Get(j);
                    if (pStr == nullptr || gStr == nullptr) {
                        return false;
                    }
                    auto itn = itsn->second.find(pStr->c_str());
                    if (itn != itsn->second.end()) {
                        if (itn->second != gStr->c_str()) {
                            result = false;
                            break;
                        }
                    } else {
                        itsn->second.emplace(pStr->c_str(), gStr->c_str());
                    }
                }
            }
            if (!result) {
                return false;
            }
        }
        return true;
    }

    void SchemeDefinitions::addFFScheme(const string& name, const map<string, string>& schemes) {
        auto it = _formFactorSchemeNames.find(name);
        if (it != _formFactorSchemeNames.end()) {
            it->second = schemes;
        } else {
            _formFactorSchemeNames.insert({name, schemes});
        }
    }

    void SchemeDefinitions::removeFFScheme(const string& name) {
        auto it = _formFactorSchemeNames.find(name);
        if (it != _formFactorSchemeNames.end()) {
            _formFactorSchemeNames.erase(it);
        }
    }

    void SchemeDefinitions::setFFInputScheme(const map<string, string>& schemes) {
        _formFactorBase = schemes;
    }

    vector<string> SchemeDefinitions::getFFSchemeNames() const {
        vector<string> names;
        for (const auto& elem : _formFactorSchemeNames) {
            names.push_back(elem.first);
        }
        return names;
    }

    map<HashId, pair<string, string>> SchemeDefinitions::getScheme(const string& name) const {
        const map<string, string>* toProcess;
        if (name.empty()) {
            toProcess = &_formFactorBase;
        } else {
            auto it = _formFactorSchemeNames.find(name);
            if (it != _formFactorSchemeNames.end()) {
                toProcess = &(it->second);
            } else {
                return map<HashId, pair<string, string>>{};
            }
        }
        map<HashId, pair<string, string>> result;
        const PID& pdg = PID::instance();
        for (const auto& elem : *toProcess) {
            vector<HashId> tmp = pdg.expandToValidVertexUIDs(elem.first, true);
            for (auto elem2 : tmp) {
                result[elem2] = make_pair(elem.first, elem.second);
            }
        }
        return result;
    }

    map<HashId, map<string, vector<string>>> SchemeDefinitions::getFFDuplicates() const {
        map<HashId, map<string, vector<string>>> result;
        const PID& pdg = PID::instance();
        for (const auto& elem : _formFactorSchemeNames) {
            for (const auto& elem2 : elem.second) {
                vector<string> chunks;
                boost::algorithm::split(chunks, elem2.second, [](char c) { return c == '_'; });
                if (chunks.size() >= 2) {
                    string token = chunks[1];
                    for (size_t idx = 2; idx < chunks.size(); ++idx) {
                        token += "_" + chunks[idx];
                    }
                    vector<HashId> tmp = pdg.expandToValidVertexUIDs(elem2.first, true);
                    for (auto id : tmp) {
                        result[id][chunks[0]].push_back(token);
                    }
                }
            }
        }
        for (const auto& elem2 : _formFactorBase) {
            vector<string> chunks;
            boost::algorithm::split(chunks, elem2.second, boost::algorithm::is_any_of("_"));
            if (chunks.size() >= 2) {
                string token = chunks[1];
                for (size_t idx = 2; idx < chunks.size(); ++idx) {
                    token += "_" + chunks[idx];
                }
                vector<HashId> tmp = pdg.expandToValidVertexUIDs(elem2.first, true);
                for (auto id : tmp) {
                    result[id][chunks[0]].push_back(token);
                }
            }
        }
        return result;
    }

    const SchemeDict<map<HashId, FFIndex>>& SchemeDefinitions::getSchemeDefs() const {
        return _formFactorSchemes;
    }

    FFIndex SchemeDefinitions::getDenominatorFormFactor(HashId processId) const {
        const auto& den = _formFactorSchemes.find("Denominator")->second;
        auto it = den.find(processId);
        if (it != den.end()) {
            return it->second;
        }
        return 0;
    }

    set<FFIndex> SchemeDefinitions::getFormFactorIndices(HashId processId) const {
        set<FFIndex> result;
        for (const auto& elem : _formFactorSchemes) {
            auto it = elem.second.find(processId);
            if (it != elem.second.end()) {
                result.insert(it->second);
            }
        }
        return result;
    }

    SchemeDict<FFIndex> SchemeDefinitions::getFFSchemesForProcess(HashId id) const {
        SchemeDict<FFIndex> result;
        for (const auto& elem : _formFactorSchemes) {
            auto it = elem.second.find(id);
            if (it != elem.second.end()) {
                result.insert({elem.first, it->second});
            } else {
                MSG_ERROR("Process not found for scheme '" + elem.first + "', hash Id: " + to_string(id));
            }
        }
        return result;
    }

    void SchemeDefinitions::init(const map<HashId, vector<string>>& formFactGroups) {
        _formFactorSchemes.clear();
        auto names = getFFSchemeNames();
        for (auto& elem : names) {
            auto res = _formFactorSchemes.insert({elem, map<HashId, FFIndex>{}});
            if (res.second) {
                map<pair<string, string>, bool> isKnownParam{};
                auto dict = getScheme(elem); // map<HashId, pair<transition, param>>
                // map<HashId, string>::iterator ite;
                for (auto& elem2 : dict) {
                    auto itp = isKnownParam.find(elem2.second);
                    if (itp == isKnownParam.end()) {
                        isKnownParam.insert({elem2.second, false});
                    }
                    auto it = formFactGroups.find(elem2.first);
                    if (it != formFactGroups.end()) {
                        auto it2 = find(it->second.begin(), it->second.end(), elem2.second.second);
                        //                        if (it2 == it->second.end()) {
                        //                            MSG_ERROR("The parametrization '" + elem2.second + "' in scheme '"
                        //                            + elem +
                        //                                      "' is unknown or misassigned. You're gonna need a bigger
                        //                                      boat.");
                        //                        }
                        if (it2 != it->second.end()) {
                            isKnownParam[elem2.second] = true;
                        }
                        ptrdiff_t pos = distance(it->second.begin(),
                                                 find(it->second.begin(), it->second.end(), elem2.second.second));
                        auto dim = static_cast<ptrdiff_t>(it->second.size());
                        if (pos >= 0 && pos < dim) {
                            res.first->second.insert({elem2.first, static_cast<size_t>(pos)});
                        }
                    }
                }
                for (auto& param : isKnownParam) {
                    if (!param.second && !param.first.second.empty()) { // allow null params for pure PS
                        MSG_ERROR("The parametrization '" + param.first.second + "' for transition '" +
                                  param.first.first + "' in scheme '" + elem +
                                  "' is unknown or misassigned. You're gonna need a bigger boat.");
                    }
                }
            }
        }
        map<pair<string, string>, bool> isKnownParam{};
        auto dict = getScheme();
        if (!dict.empty()) {
            auto res = _formFactorSchemes.insert({"Denominator", map<HashId, FFIndex>{}});
            if (res.second) {
                for (auto& elem2 : dict) {
                    auto itp = isKnownParam.find(elem2.second);
                    if (itp == isKnownParam.end()) {
                        isKnownParam.insert({elem2.second, false});
                    }
                    auto it = formFactGroups.find(elem2.first);
                    if (it != formFactGroups.end()) {
                        auto it2 = find(it->second.begin(), it->second.end(), elem2.second.second);
                        if (it2 != it->second.end()) {
                            isKnownParam[elem2.second] = true;
                        }
                        ptrdiff_t pos = distance(it->second.begin(),
                                                 find(it->second.begin(), it->second.end(), elem2.second.second));
                        auto dim = static_cast<ptrdiff_t>(it->second.size());
                        if (pos >= 0 && pos < dim) {
                            res.first->second.insert({elem2.first, static_cast<size_t>(pos)});
                        }
                    }
                }
            }
        }
        for (auto& param : isKnownParam) {
            if (!param.second && !param.first.second.empty()) { // allow null params for pure PS){
                MSG_ERROR("The parametrization '" + param.first.second + "' for transition '" + param.first.first +
                          "' in Denominator scheme is unknown or misassigned. You're gonna need a bigger boat.");
            }
        }
    }

    YAML::Emitter& operator<<(YAML::Emitter& out, const SchemeDefinitions& s) {
        out << YAML::convert<SchemeDefinitions>::encode(s);
        return out;
    }

} // namespace Hammer

namespace YAML {

    Node convert<::Hammer::SchemeDefinitions>::encode(const ::Hammer::SchemeDefinitions& value) {
        YAML::Node node;
        if (!value._formFactorSchemeNames.empty()) {
            YAML::Node numNode;
            for (const auto& elem : value._formFactorSchemeNames) {
                YAML::Node tmpAssoc;
                for (const auto& elem2 : elem.second) {
                    tmpAssoc[elem2.first] = elem2.second;
                }
                numNode[elem.first] = tmpAssoc;
            }
            node["NumeratorSchemes"] = numNode;
        }
        YAML::Node tmpAssoc;
        for (const auto& elem : value._formFactorBase) {
            tmpAssoc[elem.first] = elem.second;
        }
        node["Denominator"] = tmpAssoc;
        return node;
    }

    bool convert<::Hammer::SchemeDefinitions>::decode(const Node& node, ::Hammer::SchemeDefinitions& value) {
        if (node.IsMap()) {
            for (const auto& entry2 : node) {
                auto what = entry2.first.as<string>();
                if (what == "NumeratorSchemes") {
                    YAML::Node schemes = entry2.second;
                    if (schemes.IsMap()) {
                        value._formFactorSchemeNames.clear();
                        for (const auto& scheme : schemes) {
                            auto name = scheme.first.as<string>();
                            auto res = value._formFactorSchemeNames.insert({name, map<string, string>{}});
                            if (scheme.second.IsMap()) {
                                YAML::Node processes = scheme.second;
                                for (const auto& process : processes) {
                                    auto procname = process.first.as<string>();
                                    auto proctype = process.second.as<string>();
                                    res.first->second.insert({procname, proctype});
                                }
                            } else {
                                return false;
                            }
                        }
                    } else {
                        return false;
                    }
                } else if (what == "Denominator") {
                    value._formFactorBase.clear();
                    if (entry2.second.IsMap()) {
                        YAML::Node processes = entry2.second;
                        for (const auto& process : processes) {
                            auto procname = process.first.as<string>();
                            auto proctype = process.second.as<string>();
                            value._formFactorBase.insert({procname, proctype});
                        }
                    } else {
                        return false;
                    }
                } else {
                    return false;
                }
            }
        } else {
            return false;
        }
        return true;
    }

} // namespace YAML
