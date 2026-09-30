///
/// @file  SchemeDefinitions.hh
/// @brief Container class for Scheme Definitions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_SCHEMEDEFINITIONS_HH
#define HAMMER_SCHEMEDEFINITIONS_HH

#ifndef SCHEMEDEFINITIONS_FRIENDS
#define SCHEMEDEFINITIONS_FRIENDS typedef bool DummyFriend
#endif

#include <map>
#include <set>
#include <string>
#include <utility>
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
    class SchemeDefinitions {

        SCHEMEDEFINITIONS_FRIENDS;

    public:

        /// @name Constructors
        //@{

        /// default constructor
        SchemeDefinitions() = default;
        SchemeDefinitions(const SchemeDefinitions& other) = delete;
        SchemeDefinitions& operator=(const SchemeDefinitions& other) = delete;
        SchemeDefinitions(SchemeDefinitions&& other) = default;
        SchemeDefinitions& operator=(SchemeDefinitions&& other) = default;
        virtual ~SchemeDefinitions() noexcept = default;
        //@}

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    public:

        void write(flatbuffers::FlatBufferBuilder* msgwriter,
                   std::vector<flatbuffers::Offset<Serial::FBFFScheme>>* schemes) const;

        bool read(const Serial::FBHeader* msgreader, bool merge);

        /// @brief
        /// @param[in] schemeName
        /// @param[in] schemes
        void addFFScheme(const std::string& schemeName, const std::map<std::string, std::string>& schemes);

        /// @brief
        /// @param[in] schemeName
        void removeFFScheme(const std::string& schemeName);

        /// @brief
        /// @return
        [[nodiscard]] SchemeNameList getFFSchemeNames() const;

        /// @brief
        /// @param[in] schemes
        void setFFInputScheme(const std::map<std::string, std::string>& schemes);

        using NewLabelsByGroup = std::map<std::string, std::vector<std::string>>;

        /// @brief
        /// @return
        [[nodiscard]] HadronicIdDict<NewLabelsByGroup> getFFDuplicates() const;

    protected:

        /// @brief
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] HadronicIdDict<std::pair<std::string, std::string>>
        getScheme(const std::string& schemeName = "") const;

    public:

        /// @brief
        /// @return
        [[nodiscard]] virtual const SchemeDict<ProcIdDict<FFIndex>>& getSchemeDefs() const;

        /// @brief
        /// @param[in] id
        /// @return
        [[nodiscard]] virtual SchemeDict<FFIndex> getFFSchemesForProcess(HadronicUID id) const;

        /// @brief
        /// @param[in] processId
        /// @return
        [[nodiscard]] virtual std::set<FFIndex> getFormFactorIndices(HadronicUID processId) const;

        /// @brief
        /// @param[in] processId
        /// @return
        [[nodiscard]] virtual FFIndex getDenominatorFormFactor(HadronicUID processId) const;

        void init(const std::map<HashId, std::vector<std::string>>& formFactGroups);

    private:

        friend struct YAML::convert<::Hammer::SchemeDefinitions>;

        /// the Hammer numerator form factor schemes
        SchemeDict<std::map<std::string, std::string>> _formFactorSchemeNames;

        /// the Hammer denominator form factor scheme
        std::map<std::string, std::string> _formFactorBase;


        /// @brief Selected formfactors by subprocess
        /// @details This container provides an identifier, FFIndex, for each FF parametrization in each subprocess: A
        /// map {name of an FF scheme, map {each subprocess FF HashId (taken from _formFactors), FFIndex}}. FFIndex is
        /// the ordered position of the parametrization in the list of chosen FF parametrizations.
        SchemeDict<ProcIdDict<FFIndex>> _formFactorSchemes;
    };

    YAML::Emitter& operator<<(YAML::Emitter& out, const SchemeDefinitions& s);

} // namespace Hammer

namespace YAML {

    template <>
    struct convert<::Hammer::SchemeDefinitions> {

        static Node encode(const ::Hammer::SchemeDefinitions& value);

        static bool decode(const Node& node, ::Hammer::SchemeDefinitions& value);
    };

} // namespace YAML

#endif
