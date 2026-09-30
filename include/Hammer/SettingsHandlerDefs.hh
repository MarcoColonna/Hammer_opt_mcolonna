#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
///
/// @file  SettingsHandlerDefs.hh
/// @brief Hammer settings manager class template methods definitions
///
#pragma clang diagnostic pop

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-

#include "Hammer/Math/Utils.hh"

namespace Hammer {


    template <typename T>
    T* SettingsHandler::getNamedSettingValue(const std::string& path, const std::string& name, WTerm group) {
        auto* candidate = getEntry(path, name, group);
        if (candidate != nullptr) {
            return candidate->getValue<T>();
        }
        if (group != WTerm::COMMON) {
            return getNamedSettingValue<T>(path, name, WTerm::COMMON);
        }
        return nullptr;
    }

    template <typename T>
    T* SettingsHandler::getNamedSettingValue(const std::string& fullName, WTerm group) {
        auto pos = fullName.find(':');
        auto path = fullName.substr(0, pos);
        auto name = fullName.substr(pos + 1);
        return getNamedSettingValue<T>(path, name, group);
    }

    template <typename T>
    Setting* SettingsHandler::addSetting(const std::string& path, const std::string& name, const T& value,
                                         WTerm group) {
        auto* candidate = getEntry(path, name, group);
        if (candidate != nullptr) {
            candidate->setDefault<T>(value);
            return candidate;
        }
        Setting setitem;
        setitem.setValue<T>(value);
        setitem.setDefault<T>(value);
        return &(_settings[group][path][name] = setitem);
    }

    template <typename T>
    Setting* SettingsHandler::addSetting(const std::string& fullName, const T& value, WTerm group) {
        auto pos = fullName.find(':');
        auto path = fullName.substr(0, pos);
        auto name = fullName.substr(pos + 1);
        return addSetting<T>(path, name, value, group);
    }

    template <typename T>
    typename std::enable_if_t<std::is_same_v<typename std::remove_cv_t<T>, double>, bool>
    SettingsHandler::maybeUpdateMassWidth(const std::string& path, const std::string& name, const T& value,
                                          WTerm group) {
        if (group == WTerm::COMMON && isMassWidth(path)) {
            double tmp = toDouble(value);
            processMassWidth(path, name, tmp);
            return true;
        }
        return false;
    }

    template <typename T>
    typename std::enable_if_t<!std::is_same_v<typename std::remove_cv_t<T>, double>, bool>
    SettingsHandler::maybeUpdateMassWidth(const std::string& /*unused*/, const std::string& /*unused*/,
                                          const T& /*unused*/, WTerm /*unused*/) {
        return false;
    }

    template <typename T>
    Setting* SettingsHandler::changeSetting(const std::string& path, const std::string& name, const T& value,
                                            WTerm group) {
        if (maybeUpdateMassWidth<T>(path, name, value, group)) {
            return nullptr;
        }
        auto* candidate = getEntry(path, name, group);
        if (candidate != nullptr) {
            candidate->setValue<T>(value);
            return candidate;
        }
        if (group != WTerm::COMMON) {
            candidate = getEntry(path, name, WTerm::COMMON);
            if (candidate != nullptr) {
                auto* newset = cloneSetting(path, name, *candidate, group);
                newset->setValue<T>(value);
                return newset;
            }
        }
        return nullptr;
    }

    template <typename T>
    Setting* SettingsHandler::changeSetting(const std::string& fullName, const T& value, WTerm group) {
        auto pos = fullName.find(':');
        auto path = fullName.substr(0, pos);
        auto name = fullName.substr(pos + 1);
        return changeSetting<T>(path, name, value, group);
    }

} // namespace Hammer
