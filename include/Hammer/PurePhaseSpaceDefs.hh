///
/// @file  PurePhaseSpaceDefs.hh
/// @brief Container class for pure phase space vertices definitions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PUREPHASEPACEDEFS_HH
#define HAMMER_PUREPHASEPACEDEFS_HH

#include <set>
#include <string>
#include <vector>

#include <boost/variant.hpp>

#include "yaml-cpp/yaml.h"

#include "Hammer/IndexTypes.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/Pdg.fhh"

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
    class PurePhaseSpaceDefs {

    public:

        /// @name Constructors
        //@{

        /// default constructor
        PurePhaseSpaceDefs() = default;

        PurePhaseSpaceDefs(const PurePhaseSpaceDefs& other) = delete;
        PurePhaseSpaceDefs& operator=(const PurePhaseSpaceDefs& other) = delete;
        PurePhaseSpaceDefs(PurePhaseSpaceDefs&& other) = default;
        PurePhaseSpaceDefs& operator=(PurePhaseSpaceDefs&& other) = default;

        virtual ~PurePhaseSpaceDefs() noexcept = default;

        //@}

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    public:

        void write(flatbuffers::FlatBufferBuilder* msgwriter,
                   std::vector<flatbuffers::Offset<Serial::FBPurePS>>* purePSs) const;

        bool read(const Serial::FBHeader* msgreader, bool merge);

        /// @brief
        /// @param[in] decays
        /// @param[in] what
        void addPurePhaseSpaceVertices(const std::set<VertexName>& decays, WTerm what);

        void addPurePhaseSpaceVertex(const VertexName& decay, WTerm what);

        void clearPurePhaseSpaceVertices(WTerm what);

        /// @brief
        /// @param[in] parent
        /// @param[in] daughters
        /// @return
        [[nodiscard]] virtual NumDenPair<bool> isPurePhaseSpace(PdgId parent,
                                                                const std::vector<PdgId>& daughters) const;


    private:

        [[nodiscard]] NumDenPair<VertexUIDSet> purePhaseSpaceVertices() const;


    public:

        void init();

    private:

        friend struct YAML::convert<::Hammer::PurePhaseSpaceDefs>;

        NumDenPair<std::vector<std::string>> _purePhaseSpaceVerticesNames;

        /// @brief Pure PS vertices
        /// @details A pair of the {set of the vertex (parent + daughter) HashIds for vertices set to PS in numerator,
        /// set of the vertex (parent + daughter) HashIds for vertices set to PS in denominator}
        NumDenPair<std::set<HashId>> _purePhaseSpaceVertices;
    };

    YAML::Emitter& operator<<(YAML::Emitter& out, const PurePhaseSpaceDefs& s);

} // namespace Hammer

namespace YAML {

    template <>
    struct convert<::Hammer::PurePhaseSpaceDefs> {

        static Node encode(const ::Hammer::PurePhaseSpaceDefs& value);

        static bool decode(const Node& node, ::Hammer::PurePhaseSpaceDefs& value);
    };

} // namespace YAML


#endif
