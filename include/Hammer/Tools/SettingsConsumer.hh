///
/// @file  SettingsConsumer.hh
/// @brief Base class for accessing Hammer settings repository
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_SettingsConsumer_HH
#define HAMMER_SettingsConsumer_HH

#include <string>
#include <map>

#include "Hammer/SettingsHandler.hh"
#include "Hammer/Tools/Logging.hh"

namespace Hammer {

    /// @brief Base class to access the settings repository
    ///
    /// Implements the basic interface for querying and defining settings
    ///
    /// @ingroup Tools
    class SettingsConsumer {

    public:

        /// @name Constructors
        //@{

        /// @brief base constructor
        SettingsConsumer() = default;

        /// @brief default copy constructor
        SettingsConsumer(const SettingsConsumer&) = default;

        SettingsConsumer& operator=(const SettingsConsumer&) = default;

        /// @brief virtual destructor
        virtual ~SettingsConsumer() noexcept = default;

        //@}


        /// @brief set link to settings repository handler.
        // Static version to allow member variable initialization before any derived
        // class constructors
        /// @param[in] sh  the settings handler
        virtual void setSettingsHandler(SettingsHandler& sh);

        virtual void setSettingsHandler(const SettingsConsumer& other);

        /// @brief provide the pointer to the repository handler
        /// @return  the settings handler
        [[nodiscard]] SettingsHandler* getSettingsHandler() const;

        WTerm setWeightTerm(WTerm group);

        virtual void addRefs() const;

    protected:

        /// @brief purely virtual function for a class to define new settings
        virtual void defineSettings() = 0;


        /// @brief logging facility
        /// @return   stream to be used for logging
        [[nodiscard]] static Log& getLog();

        void initSettings();

        /// @brief provide the basic path for the settings defined by this class, as in "<path>:<setting>"
        /// @param[in] path  the class name/path
        void setPath(const std::string& path);

        /// @brief template method to add a new setting
        /// @param[in] name         the name of the setting
        /// @param[in] defaultValue default value of the setting
        template <typename T>
        void addSetting(const std::string& name, const T& defaultValue);

        /// @brief
        /// @param[in] name
        void removeSetting(const std::string& name);

        /// @brief template method to return a setting defined by this class
        /// @param[in] name         the name of the setting
        /// @return   a pointer to the setting value (null if the setting is not present or is of a different type)
        template <typename T>
        T* getSetting(const std::string& name) const;

        /// @brief template method to return a setting defined by this class
        /// @param[in] otherPath         the name of the setting
        /// @param[in] name         the name of the setting
        /// @return   a pointer to the setting value (null if the setting is not present or is of a different type)
        template <typename T>
        T* getSetting(const std::string& otherPath, const std::string& name) const;

        /// @brief method to check a boolean setting defined by this class
        /// @param[in] name         the name of the setting
        /// @return   the status of the setting
        [[nodiscard]] bool isOn(const std::string& name) const;

        /// @brief method to check a boolean setting by its name and path
        /// @param[in] otherPath    the path of the setting
        /// @param[in] name         the name of the setting
        /// @return   the status of the setting
        [[nodiscard]] bool isOn(const std::string& otherPath, const std::string& name) const;


        template <typename T>
        void updateVectorOfSettings(const std::vector<T>& values, const std::vector<std::string>& names,
                                    const std::string& path = "", WTerm group = WTerm::COMMON);

        template <typename T>
        void updateVectorOfSettings(const std::map<std::string, T>& values, const std::string& path = "",
                                    WTerm group = WTerm::COMMON);

        template <typename T>
        std::map<std::string, T> retrieveVectorOfSettings(const std::vector<std::string>& names,
                                                          const std::string& path = "",
                                                          WTerm group = WTerm::COMMON) const;

    private:

        void checkAndCleanSettings();

    protected:

        /// pointer to the settings repository handler. Static allows to be initialized once (usually in `Run` class)
        /// and shared among all other `SettingsConsumer` derived classes
        SettingsHandler* _settingHandler{nullptr};

        /// path of the settings defined by this class
        std::string _settingPath;

        WTerm _group{WTerm::COMMON};

        std::set<std::string> _ownedSettings;
    };

} // namespace Hammer

#include "Hammer/Tools/SettingsConsumerDefs.hh"

#endif
