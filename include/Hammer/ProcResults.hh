///
/// @file  ProcResults.hh
/// @brief Container for process-related results of weight calculation
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROCESSRESULTS_HH
#define HAMMER_PROCESSRESULTS_HH

#include <vector>
#include <string>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/HammerSerial.fhh"

namespace Hammer {

    class Log;
    class SpecializationDefinitions;

    /// @brief Decay process class
    ///
    /// Contains the amplitudes, weights and info associated to a decay, ...
    ///
    /// @ingroup Core
    class ProcResults : public SettingsConsumer {

    public:

        ProcResults() = default;

        explicit ProcResults(const Serial::FBProcData* msgreader);

        ProcResults(const ProcResults& other) = delete;
        ProcResults& operator=(const ProcResults& other) = delete;
        ProcResults(ProcResults&& other) = delete;
        ProcResults& operator=(ProcResults&& other) = delete;

        ~ProcResults() noexcept override = default;

        /// @brief
        /// @param[in] what
        /// @return
        const Tensor& processAmplitude(WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] what
        /// @return
        Tensor& processAmplitude(WTerm what = WTerm::NUMERATOR);

        /// @brief
        /// @param[in] what
        /// @return
        const Tensor& processAmplitudeSquared(WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] what
        /// @return
        Tensor& processAmplitudeSquared(WTerm what = WTerm::NUMERATOR);

        std::vector<std::reference_wrapper<const Tensor>> processFormFactors(const std::string& schemeName) const;

        std::vector<std::reference_wrapper<Tensor>> processFormFactors(const std::string& schemeName);

        void appendFormFactor(const std::string& schemeName, const Tensor& formfact);

        void appendFormFactor(const std::string& schemeName, Tensor&& formfact);

        void clearFormFactors();

        bool haveFormFactors() const;

        /// @brief
        /// @param[in] schemeName
        /// @return
        const Tensor& processWeight(const std::string& schemeName,
                                    const std::string& specializationId = Spec::none()) const;

        /// @brief
        /// @param[in] schemeName
        /// @return
        Tensor& processWeight(const std::string& schemeName, const std::string& specializationId = Spec::none());

        bool haveWeight(const std::string& schemeName, const std::string& specializationId = Spec::none()) const;

        void setProcessWeight(const std::string& schemeName, const Tensor& weight,
                              const std::string& specializationId = Spec::none());

        void setProcessWeight(const std::string& schemeName, Tensor&& weight,
                              const std::string& specializationId = Spec::none());

        void clearWeights();

        void removeSpecializedWeights(const std::string& specializationId);

        void removeAllSpecializedWeights();

        bool haveWeights() const;

        std::vector<std::string> availableSchemes() const;

        void setSpecDefinitions(SpecializationDefinitions* specs);

        /// @brief
        /// @param[in] msgwriter
        /// @param[in] msg
        void write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBProcData>* msg) const;

        /// @brief
        /// @param[in] msgreader
        /// @param[in] merge
        bool read(const Serial::FBProcData* msgreader, bool merge);

    protected:

        /// @brief purely virtual function for a class to define new settings
        void defineSettings() override;

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    private:


        SpecializationDict<SchemeDict<Tensor>> _processWeights;
        NumDenPair<Tensor> _processSquaredAmplitude;
        NumDenPair<Tensor> _processAmplitude;
        SchemeDict<std::vector<Tensor>> _processFormFactors;

        SpecializationDefinitions* _specializations;
    };

} // namespace Hammer

#endif
