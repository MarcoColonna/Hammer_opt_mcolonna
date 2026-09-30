///
/// @file  SpecializationDef.cc
/// @brief Hammer base amplitude class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/SpecializationDefinitions.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Math/MultiDim/ScalarContainer.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Utils.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    Log& SpecializationDefinitions::getLog() {
        return Log::getLog("Hammer.SpecializationDefinitions");
    }

    SpecializationDefinitions::SpecializationDefinitions(const IWCFFErrProviders* providers) : _providers{providers} {
    }

    void SpecializationDefinitions::write(flatbuffers::FlatBufferBuilder* msgwriter,
                                          vector<flatbuffers::Offset<Serial::FBSpecialization>>* msgs) const {
        for (const auto& elem : _definitions) {
            for (const auto& elem2 : elem.second) {
                flatbuffers::Offset<Serial::FBSpecialization> val;
                elem2.second.write(msgwriter, &val);
                msgs->push_back(val);
            }
        }
    }

    bool SpecializationDefinitions::read(const Serial::FBHeader* msgreader, bool merge) {
        if (msgreader != nullptr) {
            if (!merge) {
                _definitions.clear();
            }
            const auto* defs = msgreader->specs();
            if (defs == nullptr) {
                return false;
            }
            for (unsigned int i = 0; i < defs->size(); ++i) {
                auto item = WCSpecialization{defs->Get(i)};
                if (merge) {
                    auto it = _definitions.find(item.getPrefixId().id);
                    if (it != _definitions.end()) {
                        auto it2 = it->second.find(item.getBaseLabel());
                        if (it2 != it->second.end() && it2->second != item) {
                            MSG_ERROR(
                                "Reading mismatch specialization definition for same prefix/id combination. Cannot "
                                "merge!");
                            return false;
                        }
                    }
                }
                item.setAmplitude(_providers->getWCProvider(item.getPrefixId().prefix));
                _definitions[item.getPrefixId().id][item.getBaseLabel()] = std::move(item);
            }
            return true;
        }
        return false;
    }


    const map<IndexLabel, WCSpecialization>& SpecializationDefinitions::getSpecializations(const string& name) const {
        auto it = _definitions.find(name);
        if (it != _definitions.end()) {
            return it->second;
        }
        throw Error("Specialization name not found!");
    }

    const WCSpecialization& SpecializationDefinitions::getSpecialization(const string& name,
                                                                         IndexLabel baseLabel) const {
        auto it = _definitions.find(name);
        if (it != _definitions.end()) {
            auto it2 = it->second.find(baseLabel);
            if (it2 != it->second.end()) {
                return it2->second;
            }
            throw Error("Specialization label not found!");
        }
        throw Error("Specialization name not found!");
    }

    const WCSpecialization& SpecializationDefinitions::getSpecialization(const string& name,
                                                                         const string& prefix) const {
        auto label = _providers->getWCLabel(prefix);
        return getSpecialization(name, label);
    }


    const WCSpecialization& SpecializationDefinitions::getSpecialization(const SpecPrefixId& prefixId) const {
        auto label = _providers->getWCLabel(prefixId.prefix);
        return getSpecialization(prefixId.id, label);
    }

    WCSpecialization& SpecializationDefinitions::getSpecialization(const string& name, IndexLabel baseLabel) {
        auto it = _definitions.find(name);
        if (it != _definitions.end()) {
            auto it2 = it->second.find(baseLabel);
            if (it2 != it->second.end()) {
                return it2->second;
            }
            throw Error("Specialization label not found!");
        }
        throw Error("Specialization name not found!");
    }

    WCSpecialization& SpecializationDefinitions::getSpecialization(const string& name, const string& prefix) {
        auto label = _providers->getWCLabel(prefix);
        return getSpecialization(name, label);
    }

    WCSpecialization& SpecializationDefinitions::getSpecialization(const SpecPrefixId& prefixId) {
        auto label = _providers->getWCLabel(prefixId.prefix);
        return getSpecialization(prefixId.id, label);
    }

    void SpecializationDefinitions::addSpecialization(const string& prefix, const string& id,
                                                      const vector<string>& names) {
        IndexLabel label = _providers->getWCLabel(prefix);
        auto it = _definitions[id].find(label);
        if (it != _definitions[id].end()) {
            if (it->second.getCoordinates() != names) {
                MSG_ERROR(
                    "Cannot add specialization: specialization of the same name/prefix exists but has different "
                    "subspace definition!");
            } else {
                MSG_WARNING("Specialization " + id + " for Wilson coefficients " + prefix +
                            " already exists with the same subspace coordinate names.");
            }
        } else {
            _definitions[id][label] = WCSpecialization(prefix, label, id, names);
            _definitions[id][label].setAmplitude(_providers->getWCProvider(prefix));
            _definitions[id][label].initialize();
        }
    }

    void SpecializationDefinitions::removeSpecializations(const string& name) {
        _definitions.erase(name);
        dontUseInWeights(name);
    }

    void SpecializationDefinitions::removeSpecialization(const string& name, IndexLabel baseLabel) {
        auto it = _definitions.find(name);
        if (it != _definitions.end()) {
            it->second.erase(baseLabel);
        }
    }

    void SpecializationDefinitions::removeSpecialization(const string& name, const string& prefix) {
        auto label = _providers->getWCLabel(prefix);
        removeSpecialization(name, label);
    }

    bool SpecializationDefinitions::specializationExists(const string& name, IndexLabel baseLabel) const {
        return checkExistence(_definitions, name, baseLabel);
    }

    void SpecializationDefinitions::clear() {
        _definitions.clear();
        _useInWeights.clear();
    }

    set<string> SpecializationDefinitions::specializationIds() const {
        set<string> res;
        for (const auto& elem : _definitions) {
            res.insert(elem.first);
        }
        return res;
    }

    set<string> SpecializationDefinitions::specializationIds(IndexLabel base) const {
        set<string> res;
        for (const auto& elem : _definitions) {
            if (elem.second.find(base) != elem.second.end()) {
                res.insert(elem.first);
            }
        }
        return res;
    }

    void SpecializationDefinitions::useInWeights(const string& name) {
        if (_definitions.find(name) != _definitions.end()) {
            _useInWeights.insert(name);
        } else {
            MSG_ERROR("Unknown specialization name: " + name + ".");
        }
    }

    void SpecializationDefinitions::dontUseInWeights(const string& name) {
        _useInWeights.erase(name);
    }

    bool SpecializationDefinitions::isUsedInWeights(const string& name) const {
        return (_useInWeights.find(name) != _useInWeights.end());
    }

    const set<string>& SpecializationDefinitions::usedSpecializationsInWeights() const {
        return _useInWeights;
    }

    const Tensor& SpecializationDefinitions::getProjectionSquaredTensor(const string& name,
                                                                        IndexLabel baseLabel) const {
        auto& baselabdict = _projectionTensorCache[name];
        auto itl = baselabdict.find(baseLabel);
        if (itl != baselabdict.end()) {
            return itl->second;
        }
        if (!specializationExists(name, baseLabel)) {
            return baselabdict.emplace(baseLabel, Tensor{"SpecProj", MD::makeScalar(1.)}).first->second;
        }
        auto data = getSpecialization(name, baseLabel).getProjectionTensor();
        Tensor t = Tensor{"SpecProj", {{data, false}, {data, true}}};
        return baselabdict.emplace(baseLabel, t).first->second;
    }


    YAML::Emitter& operator<<(YAML::Emitter& out, const SpecializationDefinitions& s) {
        out << YAML::convert<SpecializationDefinitions>::encode(s);
        return out;
    }

    void operator>>(const YAML::Node& node, SpecializationDefinitions& s) {
        auto res = YAML::convert<SpecializationDefinitions>::decode(node, s);
        if (!res) {
            throw YAML::BadConversion({});
        }
    }

} // namespace Hammer

namespace YAML {

    Node convert<::Hammer::SpecializationDefinitions>::encode(const ::Hammer::SpecializationDefinitions& value) {
        Node node;
        if (!value._definitions.empty()) {
            for (const auto& elem : value._definitions) {
                Node tmp;
                for (const auto& elem2 : elem.second) {
                    tmp[elem2.second.getPrefixId().prefix] = elem2.second;
                }
                node[elem.first] = tmp;
            }
        }
        return node;
    }

    bool convert<::Hammer::SpecializationDefinitions>::decode(const Node& node,
                                                              ::Hammer::SpecializationDefinitions& value) {
        if (!node.IsMap() || value._providers == nullptr) {
            return false;
        }
        for (auto entry : node) {
            auto id = entry.first.as<string>();
            if (entry.second.IsMap()) {
                for (auto entry2 : entry.second) {
                    auto prefix = entry2.first.as<string>();
                    auto baseL = value._providers->getWCLabel(prefix);
                    auto tmp = ::Hammer::WCSpecialization{prefix, baseL, id, {}};
                    auto res = convert<::Hammer::WCSpecialization>::decode(entry2.second, tmp);
                    if (res) {
                        tmp.setAmplitude(value._providers->getWCProvider(prefix));
                        value._definitions[id][baseL] = std::move(tmp);
                    } else {
                        return false;
                    }
                }
            } else {
                return false;
            }
        }
        return true;
    }

} // namespace YAML
