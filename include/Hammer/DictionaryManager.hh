///
/// @file  DictionaryManager.hh
/// @brief Global container class for amplitudes, rates, FFs, data
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_DICTIONARY_MANAGER
#define HAMMER_DICTIONARY_MANAGER

#ifndef DICTMANAGER_FRIENDS
#define DICTMANAGER_FRIENDS typedef bool DummyFriend
#endif

#include <memory>

#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"

namespace YAML {

    class Node;

} // namespace YAML

namespace Hammer {

    // class RunDefinitions;
    class ProvidersRepo;
    class ExternalData;
    class ProcRates;
    class ProcessDefinitions;
    class SchemeDefinitions;
    class SpecializationDefinitions;
    class PurePhaseSpaceDefs;
    class Log;

    /// @brief Main class
    ///
    /// Contains ...
    ///
    /// @ingroup Core
    class DictionaryManager : public SettingsConsumer {

        DICTMANAGER_FRIENDS;

    public:

        DictionaryManager();

        DictionaryManager(const DictionaryManager& other) = delete;
        DictionaryManager& operator=(const DictionaryManager& other) = delete;
        DictionaryManager(DictionaryManager&& other) = delete;
        DictionaryManager& operator=(DictionaryManager&& other) = delete;

        ~DictionaryManager() noexcept override;

        void init();

        [[nodiscard]] virtual const ProvidersRepo& providers() const;

        [[nodiscard]] virtual const ExternalData& externalData() const;
        virtual ExternalData& externalData();

        [[nodiscard]] virtual const ProcRates& rates() const;
        virtual ProcRates& rates();

        [[nodiscard]] virtual const ProcessDefinitions& processDefs() const;
        virtual ProcessDefinitions& processDefs();

        [[nodiscard]] virtual const PurePhaseSpaceDefs& purePSDefs() const;
        virtual PurePhaseSpaceDefs& purePSDefs();

        [[nodiscard]] virtual const SchemeDefinitions& schemeDefs() const;
        virtual SchemeDefinitions& schemeDefs();

        [[nodiscard]] virtual const SpecializationDefinitions& specializationDefs() const;
        virtual SpecializationDefinitions& specializationDefs();

        void setSettingsHandler(SettingsHandler& sh) override;

        void write(flatbuffers::FlatBufferBuilder* msgwriter) const;

        bool read(const Serial::FBHeader* msgreader, bool merge);

        /// @brief read Hammer settings from a file
        /// @param[in] fileName  the file name
        void readDecays(const std::string& fileName);

        /// @brief read Hammer settings from a string
        /// @param[in] yamlData  the decay options
        void parseDecays(const std::string& yamlData);

        /// @brief write current Hammer settings to a file
        /// @param[in] fileName  the file name
        void saveDecays(const std::string& fileName);

    private:

        /// @brief
        /// @param[in] config
        void processDecays(const YAML::Node& config);

        /// @brief purely virtual function for a class to define new settings
        void defineSettings() override;

        static Log& getLog();

        std::unique_ptr<SchemeDefinitions> _schemeDefs;
        std::unique_ptr<ProvidersRepo> _providers;
        std::unique_ptr<SpecializationDefinitions> _specDefs;
        std::unique_ptr<ExternalData> _external;
        std::unique_ptr<ProcRates> _rates;
        std::unique_ptr<PurePhaseSpaceDefs> _purePSDefs;
        std::unique_ptr<ProcessDefinitions> _procDefs;
    };

} // namespace Hammer

#endif
