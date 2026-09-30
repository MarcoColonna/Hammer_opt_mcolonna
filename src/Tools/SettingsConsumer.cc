///
/// @file  SettingsConsumer.cc
/// @brief Base class for accessing Hammer settings repository
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Tools/SettingsConsumer.hh"

using namespace std;

namespace Hammer {


    void SettingsConsumer::setPath(const std::string& path) {
        _settingPath = path;
    }

    bool SettingsConsumer::isOn(const std::string& name) const {
        if (_settingHandler != nullptr) {
            bool* val = _settingHandler->getNamedSettingValue<bool>(_settingPath, name, _group);
            if (val != nullptr) {
                return *val;
            }
            throw invalid_argument("Setting '" + _settingPath + ":" + name + "' not found.");
        }
        throw Error("SettingsHandler pointer not set. Handler sought for Irish Setter and English Pointer.");
    }

    void SettingsConsumer::initSettings() {
        _ownedSettings.clear();
        defineSettings();
        checkAndCleanSettings();
        // addRefs();
    }

    void SettingsConsumer::addRefs() const {
    }

    bool SettingsConsumer::isOn(const std::string& otherPath, const std::string& name) const {
        if (_settingHandler != nullptr) {
            bool* val = _settingHandler->getNamedSettingValue<bool>(otherPath, name, _group);
            if (val != nullptr) {
                return *val;
            }
            throw invalid_argument("Setting '" + otherPath + ":" + name + "' not found.");
        }
        throw Error("SettingsHandler pointer not set. Handler sought for Irish Setter and English Pointer.");
    }

    SettingsHandler* SettingsConsumer::getSettingsHandler() const {
        return _settingHandler;
    }

    void SettingsConsumer::setSettingsHandler(SettingsHandler& sh) {
        bool needInit = (_settingHandler == nullptr);
        _settingHandler = &sh;
        if (needInit) {
            initSettings();
        }
    }

    void SettingsConsumer::setSettingsHandler(const SettingsConsumer& other) {
        bool needInit = (_settingHandler == nullptr);
        _settingHandler = other._settingHandler;
        if (needInit) {
            initSettings();
        }
    }

    WTerm SettingsConsumer::setWeightTerm(WTerm group) {
        WTerm oldgroup = _group;
        _group = group;
        return oldgroup;
    }

    void SettingsConsumer::removeSetting(const std::string& name) {
        if (_settingHandler != nullptr) {
            _settingHandler->removeSetting(_settingPath, name, _group);
        } else {
            throw InitializationError("SettingsHandler pointer not set.");
        }
    }

    void SettingsConsumer::checkAndCleanSettings() {
        if (_settingHandler != nullptr) {
            auto defined = _settingHandler->getSettings(_settingPath);
            if (defined.size() > _ownedSettings.size()) {
                vector<string> result{defined.size() - _ownedSettings.size()};
                set_difference(defined.begin(), defined.end(), _ownedSettings.begin(), _ownedSettings.end(),
                               result.begin());
                for (auto& elem : result) {
                    MSG_ERROR("Setting " + elem + " was defined before the initialization of " + _settingPath +
                              " and is not a valid setting. Removing. Please check syntax for typos.");
                    removeSetting(elem);
                }
            }
        }
    }

    Log& SettingsConsumer::getLog() {
        return Log::getLog("Hammer.SettingsConsumer");
    }

} // namespace Hammer
