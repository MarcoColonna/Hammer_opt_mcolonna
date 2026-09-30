///
/// @file  SettingsHandler.cc
/// @brief Hammer settings manager class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <fstream>

#include <boost/algorithm/string.hpp>

#include "yaml-cpp/yaml.h"

#include "Hammer/SettingsHandler.hh"
#include "Hammer/Tools/Setting.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Pdg.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/HammerYaml.hh"

using namespace std;

namespace Hammer {

    void SettingsHandler::reset() {
        for (auto& elemGroup : _settings) {
            for (auto& elemPath : elemGroup.second) {
                for (auto& elemName : elemPath.second) {
                    elemName.second.reset();
                }
            }
        }
    }

    Log& SettingsHandler::getLog() {
        return Log::getLog("Hammer.SettingsHandler");
    }

    Setting* SettingsHandler::cloneSetting(const string& path, const string& name, const Setting& value, WTerm group) {
        return &(_settings[group][path][name] = value);
    }

    Setting* SettingsHandler::cloneSetting(const string& fullName, const Setting& value, WTerm group) {
        auto pos = fullName.find(':');
        auto path = fullName.substr(0, pos);
        auto name = fullName.substr(pos + 1);
        return cloneSetting(path, name, value, group);
    }

    Setting* SettingsHandler::getEntry(const string& path, const string& name, WTerm group) {
        auto itw = _settings.find(group);
        if (itw != _settings.end()) {
            auto itg = itw->second.find(path);
            if (itg != itw->second.end()) {
                auto itn = itg->second.find(name);
                if (itn != itg->second.end()) {
                    return &(itn->second);
                }
            }
        }
        return nullptr;
    }

    Setting* SettingsHandler::resetSetting(const string& path, const string& name, WTerm group) {
        auto* candidate = getEntry(path, name, group);
        if (candidate != nullptr) {
            candidate->reset();
        }
        return candidate;
    }

    Setting* SettingsHandler::resetSetting(const string& fullName, WTerm group) {
        auto pos = fullName.find(':');
        auto path = fullName.substr(0, pos);
        auto name = fullName.substr(pos + 1);
        return resetSetting(path, name, group);
    }

    void SettingsHandler::removeSetting(const string& path, const string& name, WTerm group) {
        auto itw = _settings.find(group);
        if (itw != _settings.end()) {
            auto itg = itw->second.find(path);
            if (itg != itw->second.end()) {
                auto itn = itg->second.find(name);
                if (itn != itg->second.end()) {
                    itg->second.erase(itn);
                }
            }
        }
    }

    void SettingsHandler::removeSetting(const string& fullName, WTerm group) {
        auto pos = fullName.find(':');
        auto path = fullName.substr(0, pos);
        auto name = fullName.substr(pos + 1);
        removeSetting(path, name, group);
    }

    void SettingsHandler::renameSetting(const string& path, const string& name, const string& newName, WTerm group) {
        Setting* old = getEntry(path, name, group);
        if (old != nullptr) {
            _settings[group][path][newName] = *old;
            removeSetting(path, name, group);
        }
    }

    void SettingsHandler::renameSetting(const string& fullName, const string& newName, WTerm group) {
        auto pos = fullName.find(':');
        auto path = fullName.substr(0, pos);
        auto name = fullName.substr(pos + 1);
        renameSetting(path, name, newName, group);
    }

    set<string> SettingsHandler::getSettings(const string& path, WTerm group) const {
        set<string> tmp;
        auto itw = _settings.find(group);
        if (itw != _settings.end()) {
            auto itg = itw->second.find(path);
            if (itg != itw->second.end()) {
                for (const auto& elem : itg->second) {
                    tmp.insert(elem.first);
                }
            }
        }
        return tmp;
    }

    void SettingsHandler::parseSettings(const string& data) {
        YAML::Node config = YAML::Load(data);
        processSettings(config);
    }

    void SettingsHandler::readSettings(const string& name) {
        YAML::Node config = YAML::LoadFile(name);
        processSettings(config);
    }

    bool SettingsHandler::isMassWidth(const string& name) {
        return (name == "Masses") || (name == "Widths") || boost::starts_with(name, "Masses:") ||
               boost::starts_with(name, "Widths:");
    }

    void SettingsHandler::processMassWidth(const string& path, const string& name, double value) {
        if (boost::starts_with(path, "Masses")) {
            PID::instance().setMass(name, value);
        } else if (boost::starts_with(path, "Widths")) {
            PID::instance().setWidth(name, value);
        }
    }

    void SettingsHandler::processSetting(const YAML::Node& setting, const string& path, const string& name,
                                         WTerm group) {
        if (group == WTerm::COMMON && isMassWidth(path) && setting.Type() == YAML::NodeType::Scalar) {
            try {
                auto tmp = setting.as<double>();
                processMassWidth(path, name, tmp);
            } catch (YAML::Exception&) {
                MSG_ERROR("Problem parsing mass settings input");
            }
            return;
        }
        auto* candidate = getEntry(path, name, group);
        if (candidate == nullptr) {
            candidate = addSetting(path, name, boost::blank{}, group);
        }
        try {
            candidate->update(setting.as<Setting>());
        } catch (YAML::Exception&) {
            MSG_ERROR("Problem parsing settings input");
            return;
        }
    }

    void SettingsHandler::processSettings(const YAML::Node& input, WTerm group) {
        for (YAML::const_iterator it = input.begin(); it != input.end(); ++it) {
            auto base = it->first.as<string>();
            if (it->second.Type() == YAML::NodeType::Map) {
                YAML::Node subgroup = it->second;
                if (base == "Numerator" && group != WTerm::DENOMINATOR) {
                    processSettings(subgroup, WTerm::NUMERATOR);
                } else if (base == "Denominator" && group != WTerm::NUMERATOR) {
                    processSettings(subgroup, WTerm::DENOMINATOR);
                } else {
                    for (YAML::const_iterator it2 = subgroup.begin(); it2 != subgroup.end(); ++it2) {
                        auto optname = it2->first.as<string>();
                        processSetting(it2->second, base, optname, group);
                    }
                }
            } else {
                MSG_ERROR("Problem parsing settings input");
            }
        }
    }

    map<string, map<string, const Setting*>> SettingsHandler::getEntriesByGroup(WTerm group) const {
        map<string, map<string, const Setting*>> result;
        auto itg = _settings.find(group);
        if (itg != _settings.end()) {
            for (const auto& elem : itg->second) {
                if (elem.second.empty()) {
                    continue;
                }
                auto resp = result.insert({elem.first, NamedDict<const Setting*>{}});
                for (const auto& elem2 : elem.second) {
                    (resp.first->second)[elem2.first] = &(elem2.second);
                }
            }
        }
        return result;
    }

    void SettingsHandler::saveSettings(const string& name, bool useDefault) const {
        ofstream file;
        YAML::Emitter emitter(file);
        file.open(name.c_str());
        if (file.is_open()) {
            Setting::setEncodeUseDefault(useDefault);
            emitter << YAML::Comment(version() + string(" configuration settings")) << YAML::Newline;
            emitter << YAML::BeginMap;
            writeDict2(emitter, getEntriesByGroup(WTerm::COMMON));
            auto numeratorEntries = getEntriesByGroup(WTerm::NUMERATOR);
            if (!numeratorEntries.empty()) {
                emitter << YAML::Key << "Numerator";
                emitter << YAML::Value << YAML::BeginMap;
                writeDict2(emitter, numeratorEntries);
                emitter << YAML::EndMap;
            }
            auto denominatorEntries = getEntriesByGroup(WTerm::DENOMINATOR);
            if (!denominatorEntries.empty()) {
                emitter << YAML::Key << "Denominator";
                emitter << YAML::Value << YAML::BeginMap;
                writeDict2(emitter, denominatorEntries);
                emitter << YAML::EndMap;
            }
            emitter << YAML::EndMap;
            emitter << YAML::Newline;
            file.close();
        }
    }

    string SettingsHandler::buildName(const string& path, const string& name, WTerm group) {
        return groupToPrefix(group) + path + ":" + name;
    }

    tuple<string, string, WTerm> SettingsHandler::parseName(const string& fullName) {
        auto pos = fullName.find(':');
        string path = fullName.substr(0, pos);
        string name = fullName.substr(pos + 1);
        WTerm group = WTerm::COMMON;
        pos = name.find(':');
        if (pos != string::npos) {
            group = prefixToGroup(path);
            path = name.substr(0, pos);
            name = name.substr(pos + 1);
        }
        return tuple<string, string, WTerm>{path, name, group};
    }

    void SettingsHandler::write(flatbuffers::FlatBufferBuilder* msgwriter,
                                vector<flatbuffers::Offset<Serial::FBSetting>>* settings) const {
        size_t numEntries = 0;
        for (const auto& elemGroup : _settings) {
            for (const auto& elemPath : elemGroup.second) {
                size_t tmp = elemPath.second.size();
                if (tmp > 0) {
                    numEntries += tmp;
                    settings->reserve(numEntries);
                }
                for (const auto& elemName : elemPath.second) {
                    if (elemName.second.wasChanged()) {
                        auto serialname =
                            msgwriter->CreateString(buildName(elemPath.first, elemName.first, elemGroup.first));
                        auto serialvalue = elemName.second.write(msgwriter);
                        Serial::FBSettingBuilder serialSetting{*msgwriter};
                        serialSetting.add_name(serialname);
                        serialSetting.add_value_type(serialvalue.second);
                        serialSetting.add_value(serialvalue.first);
                        auto resS = serialSetting.Finish();
                        settings->push_back(resS);
                    }
                }
            }
        }
    }

    bool SettingsHandler::read(const Serial::FBHeader* msgreader, bool merge) {
        const auto* settings = msgreader->settings();
        if (settings == nullptr) {
            return false;
        }
        string path;
        string name;
        WTerm group = WTerm::COMMON;
        for (unsigned int i = 0; i < settings->size(); ++i) {
            const auto* setting = settings->Get(i);
            if (setting == nullptr || setting->name() == nullptr) {
                return false;
            }
            tie(path, name, group) = parseName(setting->name()->c_str());
            auto setVal = Setting{setting};
            auto* candidate = getEntry(path, name, group);
            if (merge && candidate != nullptr) {
                if (!setVal.isSame(*candidate)) {
                    string fullName = setting->name()->c_str();
                    MSG_WARNING("Setting '" + fullName + "' has different values in the two runs.");
                }
            } else {
                _settings[group][path][name] = setVal;
            }
        }
        return true;
    }

    string SettingsHandler::groupToPrefix(WTerm option) {
        switch (option) {
        case WTerm::COMMON:
            return "";
        case WTerm::NUMERATOR:
            return "Numerator:";
        case WTerm::DENOMINATOR:
            return "Denominator:";
        }
        return "";
    }

    WTerm SettingsHandler::prefixToGroup(const string& option) {
        if (option.size() <= 1) {
            return WTerm::COMMON;
        }
        if (option.find("Numerator") == 0ul) {
            return WTerm::NUMERATOR;
        }
        return WTerm::DENOMINATOR;
    }


    void SettingsHandler::saveReferences(const string& filename) const {
        ofstream reffile;
        reffile.open(filename.c_str());
        for (const auto& elem : _references) {
            reffile << elem.second;
            reffile << "\n";
        }
        reffile.close();
    }

    bool SettingsHandler::checkReference(const string& bibkey) {
        return _references.find(bibkey) != _references.end();
    }

    void SettingsHandler::addReference(const string& bibkey, const string& bibtex) {
        if (_references.find(bibkey) == _references.end()) {
            string tmpval = bibtex;
            boost::replace_all(tmpval, "\", ", "\",\n ");
            _references.insert({bibkey, tmpval});
        }
    }

    void SettingsHandler::clearReferences() {
        _references.clear();
    }


} // namespace Hammer
