///
/// @file  ProcessDefinitions.cc
/// @brief Container class for storing included/forbidden process info
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <boost/algorithm/string.hpp>

#include "yaml-cpp/yaml.h"

#include "Hammer/ProcessDefinitions.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

namespace Hammer {


    Log& ProcessDefinitions::getLog() {
        return Log::getLog("Hammer.ProcessDefinitions");
    }

    void ProcessDefinitions::write(flatbuffers::FlatBufferBuilder* msgwriter,
                                   flatbuffers::Offset<Serial::FBProcDefs>* procdefs) const {
        vector<flatbuffers::Offset<Serial::FBVecString>> includes;
        vector<flatbuffers::Offset<Serial::FBVecString>> forbids;
        for (const auto& elem : _includedDecays) {
            auto serialvals = msgwriter->CreateVectorOfStrings(elem);
            Serial::FBVecStringBuilder serialBuilder{*msgwriter};
            serialBuilder.add_value(serialvals);
            includes.push_back(serialBuilder.Finish());
        }
        for (const auto& elem : _forbiddenDecays) {
            auto serialvals = msgwriter->CreateVectorOfStrings(elem);
            Serial::FBVecStringBuilder serialBuilder{*msgwriter};
            serialBuilder.add_value(serialvals);
            includes.push_back(serialBuilder.Finish());
        }
        vector<flatbuffers::Offset<Serial::FBIdSet>> includeprocsvals;
        vector<uint64_t> includeprocskeys;
        includeprocsvals.reserve(_includedProcesses.size());
        includeprocskeys.reserve(_includedProcesses.size());
        for (auto& elem : _includedProcesses) {
            includeprocskeys.push_back(elem.first);
            vector<uint64_t> ids;
            ids.reserve(elem.second.size());
            copy(elem.second.begin(), elem.second.end(), back_inserter(ids));
            auto serialid = msgwriter->CreateVector(ids);
            Serial::FBIdSetBuilder serialBuilder{*msgwriter};
            serialBuilder.add_ids(serialid);
            includeprocsvals.push_back(serialBuilder.Finish());
        }
        vector<uint64_t> forbidprocids;
        forbidprocids.reserve(_forbiddenProcesses.size());
        copy(_forbiddenProcesses.begin(), _forbiddenProcesses.end(), back_inserter(forbidprocids));
        auto serialid = msgwriter->CreateVector(forbidprocids);
        Serial::FBIdSetBuilder serialBuilder{*msgwriter};
        serialBuilder.add_ids(serialid);
        auto serialforbidprocs = serialBuilder.Finish();
        auto serialincludes = msgwriter->CreateVector(includes);
        auto serialforbids = msgwriter->CreateVector(forbids);
        auto serialincludprocsvals = msgwriter->CreateVector(includeprocsvals);
        auto serialincludprocskeys = msgwriter->CreateVector(includeprocskeys);
        Serial::FBProcDefsBuilder builder{*msgwriter};
        builder.add_includenames(serialincludes);
        builder.add_forbidnames(serialforbids);
        builder.add_includeprocsverts(serialincludprocsvals);
        builder.add_includeprocsid(serialincludprocskeys);
        builder.add_forbidprocs(serialforbidprocs);
        *procdefs = builder.Finish();
    }

    bool ProcessDefinitions::read(const Serial::FBProcDefs* msgreader, bool merge) {
        if (msgreader != nullptr) {
            const auto* includenames = msgreader->includenames();
            if (includenames == nullptr) {
                return false;
            }
            // we should not process new decays if we are merging
            _includedDecays.clear();
            if (!merge) {
                _includedDecays.reserve(includenames->size());
                for (unsigned int i = 0; i < includenames->size(); ++i) {
                    const auto* inc = includenames->Get(i);
                    if (inc == nullptr || inc->value() == nullptr) {
                        return false;
                    }
                    const auto* vecstr = inc->value();
                    vector<string> tmp;
                    tmp.reserve(vecstr->size());
                    for (unsigned int j = 0; j < vecstr->size(); ++j) {
                        const auto* sval = vecstr->Get(j);
                        if (sval == nullptr) {
                            return false;
                        }
                        tmp.emplace_back(sval->c_str());
                    }
                    _includedDecays.push_back(tmp);
                }
            }
            const auto* forbidnames = msgreader->forbidnames();
            if (forbidnames == nullptr) {
                return false;
            }
            _forbiddenDecays.clear();
            if (!merge) {
                _forbiddenDecays.reserve(forbidnames->size());
                for (unsigned int i = 0; i < forbidnames->size(); ++i) {
                    const auto* forb = forbidnames->Get(i);
                    if (forb == nullptr || forb->value() == nullptr) {
                        return false;
                    }
                    const auto* vecstr = forb->value();
                    vector<string> tmp;
                    tmp.reserve(vecstr->size());
                    for (unsigned int j = 0; j < vecstr->size(); ++j) {
                        const auto* sval = vecstr->Get(j);
                        if (sval == nullptr) {
                            return false;
                        }
                        tmp.emplace_back(sval->c_str());
                    }
                    _forbiddenDecays.push_back(tmp);
                }
            }
            const auto* includeprocsverts = msgreader->includeprocsverts();
            const auto* includeprocsid = msgreader->includeprocsid();
            if (includeprocsid == nullptr || includeprocsverts == nullptr) {
                return false;
            }
            _includedProcesses.clear();
            if (!merge) {
                for (unsigned int i = 0; i < includeprocsid->size(); ++i) {
                    set<HashId> tmp;
                    const auto* vset = includeprocsverts->Get(i);
                    if (vset == nullptr || vset->ids() == nullptr) {
                        return false;
                    }
                    const auto* idset = vset->ids();
                    for (unsigned int j = 0; j < idset->size(); ++j) {
                        tmp.insert(idset->Get(j));
                    }
                    _includedProcesses[includeprocsid->Get(i)] = tmp;
                }
            }
            _forbiddenProcesses.clear();
            if (!merge) {
                const auto* forbidProcs = msgreader->forbidprocs();
                if (forbidProcs == nullptr || forbidProcs->ids() == nullptr) {
                    return false;
                }
                const auto* forbidprocs = forbidProcs->ids();
                for (unsigned int j = 0; j < forbidprocs->size(); ++j) {
                    _forbiddenProcesses.insert(forbidprocs->Get(j));
                }
            }
            return true;
        }
        return false;
    }

    void ProcessDefinitions::addIncludedDecay(const vector<string>& decay) {
        _includedDecays.push_back(decay);
    }

    void ProcessDefinitions::addForbiddenDecay(const vector<string>& decay) {
        _forbiddenDecays.push_back(decay);
    }

    vector<set<HashId>> ProcessDefinitions::includedDecaySignatures() const {
        return decaySignatures(_includedDecays);
    }

    vector<set<HashId>> ProcessDefinitions::forbiddenDecaySignatures() const {
        return decaySignatures(_forbiddenDecays);
    }

    vector<set<HashId>> ProcessDefinitions::decaySignatures(const vector<vector<string>>& container) {
        vector<set<HashId>> result;
        const PID& pdg = PID::instance();
        for (const auto& elem : container) {
            vector<set<HashId>> accum;
            for (const auto& elem2 : elem) {
                auto temp = pdg.expandToValidVertexUIDs(elem2);
                if (temp.empty()) {
                    accum.clear();
                    break;
                }
                if (accum.empty()) {
                    for (auto id : temp) {
                        accum.push_back({id});
                    }
                } else {
                    vector<set<HashId>> tempprod(accum.size() * temp.size());
                    for (size_t i = 0; i < accum.size(); ++i) {
                        for (size_t j = 0; j < temp.size(); ++j) {
                            tempprod[i * temp.size() + j] = accum[i];
                            tempprod[i * temp.size() + j].insert(temp[j]);
                        }
                    }
                    accum.swap(tempprod);
                }
            }
            result.insert(result.end(), accum.begin(), accum.end());
        }
        return result;
    }

    bool ProcessDefinitions::isProcessIncluded(set<HashId> subamplitudes, HashId processId) const {
        if (_includedDecayCombos.empty() && _forbiddenDecayCombos.empty()) {
            return true;
        }
        if (_includedProcesses.find(processId) != _includedProcesses.end()) {
            return true;
        }
        if (_forbiddenProcesses.find(processId) != _forbiddenProcesses.end()) {
            return false;
        }
        bool included = false;
        for (const auto& elem : _includedDecayCombos) {
            if (includes(subamplitudes.begin(), subamplitudes.end(), elem.begin(), elem.end())) {
                included = true;
                break;
            }
        }
        if (included) {
            bool forbidden = false;
            for (const auto& elem : _forbiddenDecayCombos) {
                if (subamplitudes.size() == elem.size() &&
                    equal(subamplitudes.begin(), subamplitudes.end(), elem.begin())) {
                    forbidden = true;
                    break;
                }
            }
            if (!forbidden) {
                _includedProcesses.insert(make_pair(processId, subamplitudes));
                return true;
            }
        }
        _forbiddenProcesses.insert(processId);
        return false;
    }

    set<HashId> ProcessDefinitions::decayStringsToProcessIds(const vector<vector<string>>& names) const {
        set<HashId> results;
        auto subamplitudes = decaySignatures(names);
        for (auto& elem : subamplitudes) {
            for (auto& elem2 : _includedProcesses) {
                if (includes(elem2.second.begin(), elem2.second.end(), elem.begin(), elem.end())) {
                    results.insert(elem2.first);
                }
            }
        }
        return results;
    }

    void ProcessDefinitions::forbidProcess(HashId processId) const {
        _forbiddenProcesses.insert(processId);
        _includedProcesses.erase(processId);
    }

    void ProcessDefinitions::init() {
        _includedDecayCombos = includedDecaySignatures();
        _forbiddenDecayCombos = forbiddenDecaySignatures();
        _includedProcesses.clear();
        _forbiddenProcesses.clear();
    }

    YAML::Emitter& operator<<(YAML::Emitter& out, const ProcessDefinitions& s) {
        out << YAML::convert<ProcessDefinitions>::encode(s);
        return out;
    }

} // namespace Hammer

namespace YAML {

    Node convert<::Hammer::ProcessDefinitions>::encode(const ::Hammer::ProcessDefinitions& value) {
        Node node;
        if (!value._includedDecays.empty()) {
            Node tmpNode;
            for (const auto& elem : value._includedDecays) {
                if (elem.size() == 1) {
                    tmpNode.push_back(elem[0]);
                } else {
                    Node tmpNode2;
                    for (const auto& elem2 : elem) {
                        tmpNode2.push_back(elem2);
                    }
                    tmpNode.push_back(tmpNode2);
                }
            }
            node["Include"] = tmpNode;
        }
        if (!value._forbiddenDecays.empty()) {
            Node tmpNode;
            for (const auto& elem : value._forbiddenDecays) {
                if (elem.size() == 1) {
                    tmpNode.push_back(elem[0]);
                } else {
                    Node tmpNode2;
                    for (const auto& elem2 : elem) {
                        tmpNode2.push_back(elem2);
                    }
                    tmpNode.push_back(tmpNode2);
                }
            }
            node["Forbid"] = tmpNode;
        }
        return node;
    }

    bool convert<::Hammer::ProcessDefinitions>::decode(const Node& node, ::Hammer::ProcessDefinitions& value) {
        if (node.IsMap()) {
            for (auto entry : node) {
                auto base = entry.first.as<string>();
                if (base == "Include" || base == "Forbid") {
                    vector<vector<string>>* container =
                        (base == "Include") ? &value._includedDecays : &value._forbiddenDecays;
                    container->clear();
                    if (entry.second.Type() == YAML::NodeType::Sequence) {
                        YAML::Node seqgroup = entry.second;
                        for (auto elem : seqgroup) {
                            if (elem.Type() == YAML::NodeType::Scalar) {
                                container->push_back({elem.as<string>()});
                            } else if (elem.Type() == YAML::NodeType::Sequence) {
                                try {
                                    auto tmpvec = elem.as<vector<string>>();
                                    container->push_back(tmpvec);
                                } catch (YAML::Exception&) {
                                    return false;
                                }
                            } else {
                                return false;
                            }
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
