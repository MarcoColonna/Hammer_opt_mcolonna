///
/// @file  SpecializationDefinitions.hh
/// @brief Class defining partial specialization subspace of defined external parameters
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_SPECIALIZATIONDEFINITIONS_HH
#define HAMMER_SPECIALIZATIONDEFINITIONS_HH

#include <map>
#include <set>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/IndexTypes.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/WCSpecialization.hh"
#include "Hammer/Math/Tensor.hh"

namespace Hammer {

    class Log;
    class IWCFFErrProviders;
    class AmplitudeBase;


    class SpecializationDefinitions {
    public:

        explicit SpecializationDefinitions(const IWCFFErrProviders* providers = nullptr);

        SpecializationDefinitions(const SpecializationDefinitions& other) = delete;
        SpecializationDefinitions& operator=(const SpecializationDefinitions& other) = delete;
        SpecializationDefinitions(SpecializationDefinitions&& other) = default;
        SpecializationDefinitions& operator=(SpecializationDefinitions&& other) = default;

        virtual ~SpecializationDefinitions() noexcept = default;
        //@}

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    public:

        void write(flatbuffers::FlatBufferBuilder* msgwriter,
                   std::vector<flatbuffers::Offset<Serial::FBSpecialization>>* msgs) const;

        bool read(const Serial::FBHeader* msgreader, bool merge);

        const std::map<IndexLabel, WCSpecialization>& getSpecializations(const std::string& name) const;

        const WCSpecialization& getSpecialization(const std::string& name, IndexLabel baseLabel) const;

        const WCSpecialization& getSpecialization(const std::string& name, const std::string& prefix) const;

        const WCSpecialization& getSpecialization(const SpecPrefixId& prefixId) const;

        WCSpecialization& getSpecialization(const std::string& name, IndexLabel baseLabel);

        WCSpecialization& getSpecialization(const std::string& name, const std::string& prefix);

        WCSpecialization& getSpecialization(const SpecPrefixId& prefixId);

        void addSpecialization(const std::string& prefix, const std::string& id,
                               const std::vector<std::string>& names = {});

        void removeSpecializations(const std::string& name);
        void removeSpecialization(const std::string& name, IndexLabel baseLabel);
        void removeSpecialization(const std::string& name, const std::string& prefix);

        bool specializationExists(const std::string& name, IndexLabel baseLabel) const;

        void clear();

        std::set<std::string> specializationIds() const;
        std::set<std::string> specializationIds(IndexLabel base) const;

        void useInWeights(const std::string& name);
        void dontUseInWeights(const std::string& name);
        bool isUsedInWeights(const std::string& name) const;
        const std::set<std::string>& usedSpecializationsInWeights() const;

        const Tensor& getProjectionSquaredTensor(const std::string& name, IndexLabel baseLabel) const;

    private:

        friend struct ::YAML::convert<::Hammer::SpecializationDefinitions>;

        const IWCFFErrProviders* _providers;

        SpecializationDict<std::map<IndexLabel, WCSpecialization>> _definitions;

        std::set<std::string> _useInWeights;

        mutable SpecializationDict<std::map<IndexLabel, Tensor>> _projectionTensorCache;
    };

    YAML::Emitter& operator<<(YAML::Emitter& out, const SpecializationDefinitions& s);

    void operator>>(const YAML::Node& node, SpecializationDefinitions& s);

} // namespace Hammer

namespace YAML {

    template <>
    struct convert<::Hammer::SpecializationDefinitions> {

        static Node encode(const ::Hammer::SpecializationDefinitions& value);

        static bool decode(const Node& node, ::Hammer::SpecializationDefinitions& value);
    };

} // namespace YAML

#endif
