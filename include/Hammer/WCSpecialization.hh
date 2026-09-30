///
/// @file  WCSpecialization.hh
/// @brief Class defining partial specialization subspace of defined external parameters
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_WCSPECIALIZATION_HH
#define HAMMER_WCSPECIALIZATION_HH

#ifndef WCSPECIALIZATION_FRIENDS
#define WCSPECIALIZATION_FRIENDS typedef bool DummyFriend
#endif

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/IndexTypes.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"


namespace Hammer {

    class Log;
    class AmplitudeBase;

    /// @brief Hammer parameter subspace class
    ///
    /// Contains ...
    ///
    /// @ingroup Core
    class WCSpecialization : public SettingsConsumer {

        WCSPECIALIZATION_FRIENDS;

    public:

        /// @name Constructors
        //@{

        WCSpecialization();

        WCSpecialization(const std::string& prefix, IndexLabel baseLabel, const std::string& id,
                         const std::vector<std::string>& names);

        explicit WCSpecialization(const Serial::FBSpecialization* msgreader);

        WCSpecialization(const WCSpecialization& other) = delete;
        WCSpecialization& operator=(const WCSpecialization& other) = delete;
        WCSpecialization(WCSpecialization&& other) = default;
        WCSpecialization& operator=(WCSpecialization&& other) = default; // NOLINT(bugprone-exception-escape)

        ~WCSpecialization() noexcept override = default;
        //@}

        void initialize();

        [[nodiscard]] std::vector<std::complex<double>>
        getWCSpecVectorFromDict(const std::map<std::string, std::complex<double>>& subDict) const;
        [[nodiscard]] std::vector<std::complex<double>> getWCSpecVectorFromSettings() const;

        void updateWCSpecSettings(const std::vector<std::complex<double>>& values);
        void updateWCSpecSettings(const std::map<std::string, std::complex<double>>& values);

        std::vector<std::complex<double>> resetWCSpecialization();

        std::map<std::string, std::complex<double>> retrieveWCSpecialization() const;

        /// @brief
        /// @param[in] values
        /// @param[in] data
        /// @return
        void updateWCSpecTensor(std::vector<std::complex<double>> values,
                                MultiDimensional::SharedTensorData& data) const;

        /// @brief
        /// @return prefix token label
        [[nodiscard]] bool isPartialSpecialization() const;

        /// @brief
        /// @return prefix token label
        [[nodiscard]] const SpecPrefixId& getPrefixId() const;

        /// @brief
        /// @return prefix token label
        [[nodiscard]] IndexLabel getFullLabel() const;

        /// @brief
        /// @return prefix token label
        [[nodiscard]] IndexLabel getBaseLabel() const;

        /// @brief
        /// @return prefix token label
        [[nodiscard]] uint32_t getPad() const;

        /// @brief
        /// @return prefix token label
        [[nodiscard]] const std::vector<std::string>& getCoordinates() const;

        /// @brief
        /// @return
        [[nodiscard]] MultiDimensional::SharedTensorData getProjectionTensor() const;

        /// @brief
        /// @param[in] proj
        /// @return
        void setSubspaceOrigin(const std::vector<std::complex<double>>& values);

        /// @brief
        /// @param[in] settings
        void setSubspaceOrigin(const std::map<std::string, std::complex<double>>& settings);

        /// @brief
        /// @param[in] proj
        /// @return
        void setSubspaceBasis(const std::vector<std::map<std::string, std::complex<double>>>& subspace);

        /// @brief
        /// @param[in] proj
        /// @return
        void setSubspaceVector(const std::string& coord, const std::vector<std::complex<double>>& values);

        /// @brief
        /// @param[in] proj
        /// @return
        void setSubspaceVector(const std::string& coord, const std::map<std::string, std::complex<double>>& values);

        /// @brief
        /// @param[in] amp
        /// @return
        void setAmplitude(AmplitudeBase* amp);

        void write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBSpecialization>* msg) const;

    protected:

        /// @brief defines new settings for this class
        void defineSettings() override;

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    private:

        friend struct ::YAML::convert<::Hammer::WCSpecialization>;

        friend bool operator==(const WCSpecialization& /*a*/, const WCSpecialization& /*b*/);

        friend bool operator!=(const WCSpecialization& /*a*/, const WCSpecialization& /*b*/);

        void read(const Serial::FBSpecialization* msgreader);

        [[nodiscard]] size_t coordIndex(const std::string& name) const;

        // void updateProjTensorInSettings();

        SpecPrefixId _prefixId;
        uint32_t _pad;
        std::vector<std::string> _coordNames;
        IndexLabel _baseLabel, _fullLabel;
        MultiDimensional::SharedTensorData _projector;

        AmplitudeBase* _baseAmpl;
    };


    bool operator==(const WCSpecialization& a, const WCSpecialization& b);

    bool operator!=(const WCSpecialization& a, const WCSpecialization& b);

} // namespace Hammer

namespace YAML {

    template <>
    struct convert<::Hammer::WCSpecialization> {

        static Node encode(const ::Hammer::WCSpecialization& value);

        static bool decode(const Node& node, ::Hammer::WCSpecialization& value);
    };

} // namespace YAML

#endif
