///
/// @file  ProcessDefinitions.hh
/// @brief Container class for storing included/forbidden process info
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROCESSDEFINITIONS_HH
#define HAMMER_PROCESSDEFINITIONS_HH

#ifndef PROCESSDEFINITIONS_FRIENDS
#define PROCESSDEFINITIONS_FRIENDS typedef bool DummyFriend
#endif

#include <set>
#include <string>
#include <vector>

#include <boost/variant.hpp>

#include "yaml-cpp/yaml.h"

#include "Hammer/IndexTypes.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/Logging.hh"

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
    class ProcessDefinitions {

        PROCESSDEFINITIONS_FRIENDS;

    public:

        /// @name Constructors
        //@{

        /// default constructor
        ProcessDefinitions() = default;

        ProcessDefinitions(const ProcessDefinitions& other) = delete;
        ProcessDefinitions& operator=(const ProcessDefinitions& other) = delete;
        ProcessDefinitions(ProcessDefinitions&& other) = default;
        ProcessDefinitions& operator=(ProcessDefinitions&& other) = default;

        virtual ~ProcessDefinitions() noexcept = default;
        //@}

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    public:

        void write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBProcDefs>* procdefs) const;

        bool read(const Serial::FBProcDefs* msgreader, bool merge);

        /// @brief
        /// @param[in] decay
        void addIncludedDecay(const std::vector<VertexName>& decay);

        /// @brief
        /// @param[in] decay
        void addForbiddenDecay(const std::vector<VertexName>& decay);

        std::set<HashId> decayStringsToProcessIds(const std::vector<std::vector<VertexName>>& names) const;

    protected:

        /// @brief
        /// @return
        std::vector<std::set<HashId>> includedDecaySignatures() const;

        /// @brief
        /// @return
        std::vector<std::set<HashId>> forbiddenDecaySignatures() const;

        /// @brief
        /// @param[in] container
        /// @return
        static std::vector<std::set<HashId>> decaySignatures(const std::vector<std::vector<std::string>>& container);

    public:

        /// @brief
        /// @param[in] subamplitudes
        /// @param[in] processId
        /// @return
        virtual bool isProcessIncluded(std::set<HashId> subamplitudes, ProcessUID processId) const;

        /// @brief
        /// @param[in] processId
        virtual void forbidProcess(ProcessUID processId) const;

        void init();

    private:

        friend struct YAML::convert<::Hammer::ProcessDefinitions>;

        /// the Hammer included decays
        std::vector<std::vector<std::string>> _includedDecays;

        /// the Hammer forbidden decays
        std::vector<std::vector<std::string>> _forbiddenDecays;

        /// @brief Included amplitude combinations
        /// @details A list of {the set of amplitude class HashIds} composing each included process specification
        std::vector<std::set<HashId>> _includedDecayCombos;
        /// @brief Forbidden amplitude combinations
        /// @details A list of {the set of amplitude class HashIds} composing each forbidden process specification
        std::vector<std::set<HashId>> _forbiddenDecayCombos;
        /// @brief Included process HashIds
        /// @details The set of {the full HashId} for each included process, generated from the HashIds of the process'
        /// component amplitudes
        mutable ProcIdDict<std::set<HashId>> _includedProcesses;
        /// @brief Forbidden process HashIds
        /// @details The set of {the full HashId} for each forbidden process, generated from the HashIds of the process'
        /// component amplitudes
        mutable std::set<HashId> _forbiddenProcesses;
    };

    YAML::Emitter& operator<<(YAML::Emitter& out, const ProcessDefinitions& s);

} // namespace Hammer

namespace YAML {

    template <>
    struct convert<::Hammer::ProcessDefinitions> {

        static Node encode(const ::Hammer::ProcessDefinitions& value);

        static bool decode(const Node& node, ::Hammer::ProcessDefinitions& value);
    };

} // namespace YAML

#endif
