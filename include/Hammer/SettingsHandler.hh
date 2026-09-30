///
/// @file  SettingsHandler.hh
/// @brief Hammer settings manager class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_SettingsHandler_HH
#define HAMMER_SettingsHandler_HH

#include <map>
#include <set>
#include <string>
#include <vector>
#include <type_traits>

#include <boost/variant.hpp>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Tools/Setting.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/Logging.hh"

namespace YAML {

    class Node;
    class Emitter;

} // namespace YAML

namespace Hammer {

    class Log;

    /// @brief Hammer settings manager class
    ///
    /// Stores Hammer options, provides facilities for saving and reading option files,
    /// provides option query interface. It also provide a repository of other run-wide
    /// information, such as the bibliography associated to the specific run, the list
    /// of files being processed
    ///
    /// @ingroup Handlers
    class SettingsHandler {

    public:

        /// @name Constructors
        //@{

        /// default constructor
        SettingsHandler() = default;

        SettingsHandler(const SettingsHandler& other) = delete;
        SettingsHandler& operator=(const SettingsHandler& other) = delete;
        SettingsHandler(SettingsHandler&& other) = delete;
        SettingsHandler& operator=(SettingsHandler&& other) = delete;

        ~SettingsHandler() noexcept = default;

        //@}

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        [[nodiscard]] static Log& getLog();

    public:

        /// @brief reset all the contents of `SettingsHandler`
        void reset();

        /// @name Hammer options
        //@{

        /// @brief access the list of registered settings
        /// @param[in]   path      the class owning the settings. An empty string returns all the registered settings
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        /// @return a list of setting names
        [[nodiscard]] std::set<std::string> getSettings(const std::string& path = "",
                                                        WTerm group = WTerm::COMMON) const;

        /// @brief access a specific setting value
        /// @param[in]   path      the class owning the settings
        /// @param[in]   name      the name of the setting
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        /// @return a pointer to the value, `nullptr` if the setting is not found or the type does not match
        template <typename T>
        T* getNamedSettingValue(const std::string& path, const std::string& name, WTerm group = WTerm::COMMON);

        /// @brief access a specific setting value
        /// @param[in]   fullName      the full name of the setting in the "<path>:<name>" format
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        /// @return a pointer to the value, `nullptr` if the setting is not found or the type does not match
        template <typename T>
        T* getNamedSettingValue(const std::string& fullName, WTerm group = WTerm::COMMON);

        /// @brief add a setting to the store
        /// @param[in]   path      the class owning the settings
        /// @param[in]   name      the name of the setting
        /// @param[in]   value     the value (and default value) of the setting. If setting is already present only
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        /// default value is changed
        /// @return the added setting
        template <typename T>
        Setting* addSetting(const std::string& path, const std::string& name, const T& value,
                            WTerm group = WTerm::COMMON);

        /// @brief add a setting to the store
        /// @param[in]   fullName      the full name of the setting in the "<path>:<name>" format
        /// @param[in]   value     the value (and default value) of the setting
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        /// @return the added setting
        template <typename T>
        Setting* addSetting(const std::string& fullName, const T& value, WTerm group = WTerm::COMMON);

        Setting* cloneSetting(const std::string& path, const std::string& name, const Setting& value,
                              WTerm group = WTerm::COMMON);

        Setting* cloneSetting(const std::string& fullName, const Setting& value, WTerm group = WTerm::COMMON);

        Setting* resetSetting(const std::string& path, const std::string& name, WTerm group = WTerm::COMMON);

        Setting* resetSetting(const std::string& fullName, WTerm group = WTerm::COMMON);

        Setting* getEntry(const std::string& path, const std::string& name, WTerm group);

        /// @brief
        /// @param[in] path
        /// @param[in] name
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        void removeSetting(const std::string& path, const std::string& name, WTerm group = WTerm::COMMON);

        /// @brief
        /// @param[in] fullName
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        void removeSetting(const std::string& fullName, WTerm group = WTerm::COMMON);

        /// @brief change a setting value
        /// @param[in]   path      the class owning the settings
        /// @param[in]   name      the name of the setting
        /// @param[in]   value     the new value of the setting
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        /// @return the changed setting, `nullptr` if the setting was not found
        template <typename T>
        Setting* changeSetting(const std::string& path, const std::string& name, const T& value,
                               WTerm group = WTerm::COMMON);

        /// @brief change a setting value
        /// @param[in]   fullName      the full name of the setting in the "<path>:<name>" format
        /// @param[in]   value     the new value of the setting
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        /// @return the changed setting, `nullptr` if the setting was not found
        template <typename T>
        Setting* changeSetting(const std::string& fullName, const T& value, WTerm group = WTerm::COMMON);

        void renameSetting(const std::string& path, const std::string& name, const std::string& newName,
                           WTerm group = WTerm::COMMON);

        void renameSetting(const std::string& fullName, const std::string& newName, WTerm group = WTerm::COMMON);

        /// @brief read Hammer settings from a file
        /// @param[in] fileName  the file name
        void readSettings(const std::string& fileName);

        /// @brief read Hammer settings from a string
        /// @param[in] data  the settings
        void parseSettings(const std::string& data);

        /// @brief write current Hammer settings to a file
        /// @param[in] name  the file name
        /// @param[in] useDefault  the file name
        void saveSettings(const std::string& name, bool useDefault = true) const;

        //@}

        void write(flatbuffers::FlatBufferBuilder* msgwriter,
                   std::vector<flatbuffers::Offset<Serial::FBSetting>>* settings) const;

        bool read(const Serial::FBHeader* msgreader, bool merge);


        /// @name Bibliography
        //@{

        /// writes a BibTeX file containing the bibliography associated to current run
        /// @param[in] filename  the name of the file
        void saveReferences(const std::string& filename) const;

        /// check a bibliographic reference to the current run
        /// @param[in] bibkey   the BibTeX key
        bool checkReference(const std::string& bibkey);

        /// add a bibliographic reference to the current run
        /// @param[in] bibkey   the BibTeX key
        /// @param[in] bibtex   the body of the citation in BibTeX format
        void addReference(const std::string& bibkey, const std::string& bibtex);

        /// empties current run bibliography
        void clearReferences();

        //@}

    private:

        /// @brief
        /// @param[in] setting
        /// @param[in] path
        /// @param[in] name
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        void processSetting(const YAML::Node& setting, const std::string& path, const std::string& name,
                            WTerm group = WTerm::COMMON);

        /// @brief
        /// @param[in] input
        /// @param[in]   group     whether it applies to the numerator, denominator or both
        void processSettings(const YAML::Node& input, WTerm group = WTerm::COMMON);

        [[nodiscard]] static std::string groupToPrefix(WTerm option);

        [[nodiscard]] static WTerm prefixToGroup(const std::string& option);

        [[nodiscard]] static bool isMassWidth(const std::string& name);
        static void processMassWidth(const std::string& path, const std::string& name, double value);

        [[nodiscard]] std::map<std::string, std::map<std::string, const Setting*>> getEntriesByGroup(WTerm group) const;

        [[nodiscard]] static std::string buildName(const std::string& path, const std::string& name, WTerm group);

        [[nodiscard]] static std::tuple<std::string, std::string, WTerm> parseName(const std::string& fullName);

        template <typename T>
        typename std::enable_if_t<std::is_same_v<typename std::remove_cv_t<T>, double>, bool>
        maybeUpdateMassWidth(const std::string& path, const std::string& name, const T& value,
                             WTerm group = WTerm::COMMON);

        template <typename T>
        typename std::enable_if_t<!std::is_same_v<typename std::remove_cv_t<T>, double>, bool>
        maybeUpdateMassWidth(const std::string& path, const std::string& name, const T& value,
                             WTerm group = WTerm::COMMON);


        template <typename T>
        using NamedDict = std::map<std::string, T>;

        /// the Hammer options dictionary
        std::map<WTerm, NamedDict<NamedDict<Setting>>> _settings;

        std::map<std::string, std::string> _references;
    };

} // namespace Hammer

#include "Hammer/SettingsHandlerDefs.hh"

#endif
