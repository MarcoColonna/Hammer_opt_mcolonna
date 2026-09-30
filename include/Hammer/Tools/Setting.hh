///
/// @file  Setting.hh
/// @brief Hammer setting class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_Setting_HH
#define HAMMER_Setting_HH

#include <complex>
#include <string>
#include <utility>
#include <vector>

#include <boost/variant.hpp>

#include "yaml-cpp/yaml.h"

#include "Hammer/Tools/HammerSerial.fhh"

namespace Hammer {

    class Log;
    struct SettingYamlConverter;
    /// @brief container for an Hammer run option
    ///
    /// stores the value, name, owner, description and default value of an Hammer option
    ///
    /// @ingroup Handlers
    class Setting {

    public:

        /// default constructor
        Setting() = default;
        Setting(const Setting& other) = default;
        Setting& operator=(const Setting& other);
        Setting(Setting&& other) = default;
        Setting& operator=(Setting&& other) = default;
        ~Setting() = default;

        explicit Setting(const Serial::FBSetting* msgreader);

        /// @brief access the setting value
        /// @return a pointer to the value, `nullptr` if the data type does not match
        template <typename T>
        T* getValue();

        /// @brief access the setting value (`const` version)
        /// @return a pointer to the value, `nullptr` if the data type does not match
        template <typename T>
        const T* getValue() const;

        /// @brief modify the setting value
        /// @param[in] value  the new value of the setting
        template <typename T>
        void setValue(const T& value);

        /// @brief get the default value of this setting
        /// @return a pointer to the default value, `nullptr` if the data type does not match
        template <typename T>
        T* getDefault();

        /// @brief get the default value of this setting (`const` version)
        /// @return a pointer to the default value, `nullptr` if the data type does not match
        template <typename T>
        const T* getDefault() const;

        /// @brief modify the setting default value
        /// @param[in] value  the new default value of the setting
        template <typename T>
        void setDefault(const T& value);

        /// @brief reset the contents of the settings
        void reset();

        [[nodiscard]] bool wasChanged() const;

        [[nodiscard]] bool isSame(const Setting& other) const;

        void setDefault();

        void update(const Setting& other);

        /// @brief properly format the value of the setting (so that can be e.g. displayed on a screen)
        /// @param[in] useDefault    whether to use the current value or the default value
        /// @return the formatted string
        [[nodiscard]] std::string toString(bool useDefault = false) const;

        using WrittenSettingType = std::pair<flatbuffers::Offset<void>, Serial::FBSettingTypes>;

        WrittenSettingType write(flatbuffers::FlatBufferBuilder* msgwriter) const;

        static void setEncodeUseDefault(bool useDefault);

    private:

        friend struct SettingEncoder;
        friend struct SettingChecker;
        friend class SettingWriter;
        friend struct SettingStringConverter;
        friend struct YAML::convert<::Hammer::Setting>;

        using MatrixType = std::vector<std::vector<double>>;

        /// data type for the setting values
        using SettingType = boost::variant<boost::blank, bool, int, double, std::string, std::complex<double>,
                                           std::vector<std::string>, std::vector<double>, MatrixType>;

        void read(const Serial::FBSetting* msgreader);


        /// setting value
        SettingType _value;
        /// setting default value
        SettingType _default;

        static bool _encodeUseDefault;
    };

    YAML::Emitter& operator<<(YAML::Emitter& out, const Setting& s);

} // namespace Hammer

namespace YAML {

    template <>
    struct convert<::Hammer::Setting> {

        static Node encode(const ::Hammer::Setting& value);

        static bool decode(const Node& node, ::Hammer::Setting& value);
    };

} // namespace YAML

#include "Hammer/Tools/SettingDefs.hh"

#endif
