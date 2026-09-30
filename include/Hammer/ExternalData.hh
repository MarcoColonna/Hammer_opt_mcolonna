///
/// @file  ExternalData.hh
/// @brief Container class for values of WC and FF vectors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_EXTERNALDATA
#define HAMMER_EXTERNALDATA

#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include <boost/thread/tss.hpp>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/Pdg.fhh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/Utils.hh"

namespace Hammer {

    class AmplitudeBase;
    class FormFactorBase;
    class Tensor;
    class RunDefinitions;
    class IWCFFErrProviders;
    class SpecializationDefinitions;

    namespace MD = MultiDimensional;

    /// @brief Main class
    ///
    /// Contains ...
    ///
    /// @ingroup Core
    class ExternalData : public SettingsConsumer {

    public:

        explicit ExternalData(const IWCFFErrProviders* provs, const SpecializationDefinitions* specs);

        ExternalData(const ExternalData& other) = delete;
        ExternalData& operator=(const ExternalData& other) = delete;
        ExternalData(ExternalData&& other) = delete;
        ExternalData& operator=(ExternalData&& other) = delete;

        ~ExternalData() noexcept override;

        /// @brief
        /// @param[in] schemeName
        /// @param[in] labels
        /// @return
        virtual const Tensor& getExternalVectors(const std::string& schemeName, LabelsList labels,
                                                 const std::string& specializationId = Spec::none()) const;

        /// @brief
        /// @param[in] parent
        /// @param[in] daughters
        /// @param[in] granddaughters
        /// @param[in] what
        /// @return
        virtual MD::SharedTensorData getWilsonCoefficients(PdgId parent, const std::vector<PdgId>& daughters,
                                                           const std::vector<PdgId>& granddaughters = {},
                                                           WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] amp
        /// @param[in] what
        /// @return
        virtual MD::SharedTensorData getWilsonCoefficients(const AmplitudeBase* amp,
                                                           WTerm what = WTerm::NUMERATOR) const;


        /// @brief
        /// @param[in] prefixName
        /// @param[in] values
        /// @param[in] what
        void setWilsonCoefficients(const std::string& prefixName, const std::vector<std::complex<double>>& values,
                                   WTerm what = WTerm::NUMERATOR, bool updateSettings = true);

        void setWilsonCoefficients(const std::string& prefixName,
                                   const std::map<std::string, std::complex<double>>& values,
                                   WTerm what = WTerm::NUMERATOR, bool updateSettings = true);

        void resetWilsonCoefficients(const std::string& prefixName, WTerm what = WTerm::NUMERATOR);

        [[nodiscard]] std::map<std::string, std::complex<double>>
        retrieveWilsonCoefficients(const std::string& prefixName, WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] ff
        /// @param[in] schemeName
        /// @return
        virtual MD::SharedTensorData getFFEigenVectors(FormFactorBase* ff, const std::string& schemeName) const;


        /// @brief
        /// @param[in] process
        /// @param[in] values
        void setFFEigenVectors(const FFPrefixGroup& process, const std::vector<double>& values,
                               bool updateSettings = true);

        /// @brief
        /// @param[in] process
        /// @param[in] values
        void setFFEigenVectors(const FFPrefixGroup& process, const std::map<std::string, double>& values,
                               bool updateSettings = true);

        void resetFFEigenVectors(const FFPrefixGroup& process);

        [[nodiscard]] std::map<std::string, double> retrieveFFEigenvectors(const FFPrefixGroup& process) const;

        MD::SharedTensorData getTempFFEigenVectors(const FFPrefixGroup& process,
                                                   const std::vector<double>& values) const;
        MD::SharedTensorData getTempFFEigenVectors(const FFPrefixGroup& process,
                                                   const std::map<std::string, double>& values) const;


        void init(std::vector<std::string> schemeNames);

        void reInitFormFactorErrors();

    private:

        void initESSettings(const SpecPrefixId& prefixid, const std::vector<std::string>& names) const;

        void initWilsonCoefficients();
        void initFormFactorErrors();
        void initExternalVectors();

        void defineSettings() override;

        static Log& getLog();

        std::vector<std::string> _schemes;

        const IWCFFErrProviders* _providers;
        const SpecializationDefinitions* _specializations;


        // all the WC vectors (in {numerator, denominator} form) that can be used to contract tensor weights.
        // They are indexed by the corresponding (positive) IndexLabel. For partial specializations,
        // the new remaining WC subspace vectors are stored here, indexed by the corresponding extended
        // IndexLabel. The usage of thread specific pointer allows having independent dictionaries per-thread,
        // (used in the *Local version)
        using processWilsonCoefficientsType = std::map<IndexLabel, std::array<MD::SharedTensorData, 2>>;
        mutable boost::thread_specific_ptr<processWilsonCoefficientsType> _processWilsonCoefficients;

        // Same, but for Form Factor EigenVectors. For each label, they are organized by scheme.
        // Only one is stored as denominators don't have eigenvectors (no need from EVTGEN etc).
        using processFFEigenVectorsType = std::map<IndexLabel, SchemeDict<MD::SharedTensorData>>;
        mutable boost::thread_specific_ptr<processFFEigenVectorsType> _processFFEigenVectors;


        // Ready to use external data organized by scheme and by IndexLabels built from the contents in
        // _processWilsonCoefficients and _processFFEigenVectors and neatly packaged into Tensor object ready
        // to contract. The list of labels correspond to all the indices that will be contracted by that Tesnor.
        // Since_processWilsonCoefficients and _processFFEigenVectors use SharedTensorData and the Tensors in
        // _externalVectors are OuterTensors, no unnecessary copies are created. Moreover changes in the underlying
        // vectors automatically propagate
        using externalVectorsType = std::unordered_map<std::string, SpecializationDict<UMap<LabelsList, Tensor>>>;
        mutable boost::thread_specific_ptr<externalVectorsType> _externalVectors;
    };

} // namespace Hammer

#endif
