#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
///
/// @file  SettingDefs.hh
/// @brief Settings class template methods definitions
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
    [[nodiscard]] T* Setting::getValue() {
        return boost::get<T>(&_value);
    }

    template <typename T>
    [[nodiscard]] const T* Setting::getValue() const {
        return boost::get<const T>(&_value);
    }

    template <typename T>
    void Setting::setValue(const T& value) {
        if (_value.which() != 0 && _value.type() != typeid(T)) {
            throw InitializationError("Invalid type in setting!");
        }
        _value = value;
    }

    template <typename T>
    T* Setting::getDefault() {
        return boost::get<T>(&_default);
    }

    template <typename T>
    [[nodiscard]] const T* Setting::getDefault() const {
        return boost::get<const T>(&_default);
    }

    template <typename T>
    void Setting::setDefault(const T& value) {
        if (_default.which() != 0 && _default.type() != typeid(T)) {
            if (_default.which() == 7 and typeid(T) == typeid(std::complex<double>)) {
                const std::vector<double>* tmp = boost::get<std::vector<double>>(&_value);
                if (tmp != nullptr && tmp->size() == 2) {
                    _default = value;
                    _value = std::complex<double>((*tmp)[0], (*tmp)[1]);
                    return;
                }
            }
            throw InitializationError("Invalid type in setting!");
        }
        _default = value;
    }

} // namespace Hammer
