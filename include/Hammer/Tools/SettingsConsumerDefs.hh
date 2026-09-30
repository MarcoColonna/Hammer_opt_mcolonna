#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
///
/// @file  SettingsConsumerDefs.hh
/// @brief Settings Consumer class template methods definitions
///
#pragma clang diagnostic pop

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include "Hammer/Exceptions.hh"

namespace Hammer {

    template <typename T>
    void SettingsConsumer::addSetting(const std::string& name, const T& defaultValue) {
        if (_settingHandler) {
            _ownedSettings.insert(name);
            _settingHandler->addSetting<T>(_settingPath, name, defaultValue, _group);
        } else {
            throw InitializationError("SettingsHandler pointer not set.");
        }
    }

    template <typename T>
    T* SettingsConsumer::getSetting(const std::string& name) const {
        if (_settingHandler) {
            return _settingHandler->getNamedSettingValue<T>(_settingPath, name, _group);
        }
        throw InitializationError("SettingsHandler pointer not set.");
    }

    template <typename T>
    T* SettingsConsumer::getSetting(const std::string& otherPath, const std::string& name) const {
        if (_settingHandler) {
            return _settingHandler->getNamedSettingValue<T>(otherPath, name, _group);
        }
        throw InitializationError("SettingsHandler pointer not set.");
    }

    template <typename T>
    void SettingsConsumer::updateVectorOfSettings(const std::vector<T>& values, const std::vector<std::string>& names,
                                                  const std::string& path, WTerm group) {
        if (_settingHandler) {
            ASSERT(values.size() == names.size());
            auto tmpPath = (!path.empty()) ? path : _settingPath;
            for (size_t pos = 0ul; pos < names.size(); ++pos) {
                _settingHandler->changeSetting<T>(tmpPath, names[pos], values[pos], group);
            }

        } else {
            throw InitializationError("SettingsHandler pointer not set.");
        }
    }

    template <typename T>
    void SettingsConsumer::updateVectorOfSettings(const std::map<std::string, T>& values, const std::string& path,
                                                  WTerm group) {
        if (_settingHandler) {
            auto tmpPath = (!path.empty()) ? path : _settingPath;
            for (auto& elem : values) {
                if (_settingHandler->changeSetting<T>(tmpPath, elem.first, elem.second, group) == nullptr) {
                    MSG_WARNING("Setting '" + elem.first + "' in path '" + path + "' does not exist.");
                }
            }
        } else {
            throw InitializationError("SettingsHandler pointer not set.");
        }
    }

    template <typename T>
    std::map<std::string, T> SettingsConsumer::retrieveVectorOfSettings(const std::vector<std::string>& names,
                                                                        const std::string& path, WTerm group) const {
        if (_settingHandler) {
            std::map<std::string, T> result;
            auto tmpPath = (!path.empty()) ? path : _settingPath;
            for (const auto& name : names) {
                auto* pValue = _settingHandler->getNamedSettingValue<T>(path, name, group);
                if (pValue != nullptr) {
                    result.emplace(name, *pValue);
                } else {
                    MSG_WARNING("Setting '" + name + "' in path '" + path + "' does not exist.");
                }
            }
            return result;
        }
        throw InitializationError("SettingsHandler pointer not set.");
    }


} // namespace Hammer
