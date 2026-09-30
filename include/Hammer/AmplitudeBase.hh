///
/// @file  AmplitudeBase.hh
/// @brief Hammer base amplitude class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPLITUDE_BASE
#define HAMMER_AMPLITUDE_BASE

#include <string>
#include <vector>

#include "Hammer/Particle.fhh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/ParticleData.hh"
#include "Hammer/Tools/SettingsConsumer.hh"

namespace Hammer {

    class Log;

    /// @brief Base class for amplitudes
    ///
    /// Provides the amplitude contents as tensor...
    ///
    /// @ingroup Base
    class AmplitudeBase : public ParticleData, public SettingsConsumer {

    public:

        AmplitudeBase() = default;

        AmplitudeBase(const AmplitudeBase& other) = delete;
        AmplitudeBase& operator=(const AmplitudeBase& other) = delete;

        AmplitudeBase(AmplitudeBase&& other) = delete;
        AmplitudeBase& operator=(AmplitudeBase&& other) = delete;
        ~AmplitudeBase() noexcept override = default;

        /// @brief method to evaluate the object on a specific particle set
        /// @param[in] parent  the parent Particle
        /// @param[in] daughters the daughters (and grand-daughters, if necessary) Particle list
        /// @param[in] references the parent Particle siblings (necessary e.g. for helicity amplitude phase conventions)
        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) override = 0;

        [[nodiscard]] std::vector<std::complex<double>>
        getWCVectorFromDict(const std::map<std::string, std::complex<double>>& wcDict) const;
        [[nodiscard]] std::map<std::string, std::complex<double>>
        getDictFromWCVector(const std::vector<std::complex<double>>& wcVect) const;
        [[nodiscard]] std::vector<std::complex<double>> getWCVectorFromSettings(WTerm what) const;
        std::vector<std::complex<double>> resetWCSettings(WTerm what);

        void updateWCSettings(const std::vector<std::complex<double>>& values, WTerm what);
        void updateWCSettings(const std::map<std::string, std::complex<double>>& values, WTerm what);

        [[nodiscard]] std::map<std::string, std::complex<double>> retrieveWCSettings(WTerm what) const;

        /// @brief
        /// @param[in] values
        /// @param[in] data
        /// @return
        void updateWCTensor(std::vector<std::complex<double>> values, MultiDimensional::SharedTensorData& data) const;

        /// @brief
        /// @return
        [[nodiscard]] std::pair<std::string, IndexLabel> getWCInfo() const;

        /// @brief generate the projection tensor for WC subspace specialization
        /// @param[in] vector of map wc_i -> a_i
        void createWCProjectionTensor(const std::vector<std::map<std::string, std::complex<double>>>& subspace,
                                      IndexLabel label, const MultiDimensional::SharedTensorData& origin,
                                      MultiDimensional::SharedTensorData& proj) const;

        void updateWCProjectionVector(int index, const std::map<std::string, std::complex<double>>& direction,
                                      MultiDimensional::SharedTensorData& proj) const;

        void updateWCProjectionVector(int index, const std::vector<std::complex<double>>& direction,
                                      MultiDimensional::SharedTensorData& proj) const;

        /// @brief initializes the amplitude (defines settings associated to this amplitude, etc.)
        void init();

        /// @brief returns a reference to itself as a Tensor
        /// @return itself
        Tensor& getTensor();

        /// @brief returns a reference to itself as a Tensor
        /// @return itself
        [[nodiscard]] const Tensor& getTensor() const;

        /// @brief select a specific signature to be the current signature
        /// @param[in] idx  the signature index
        bool setSignatureIndex(size_t idx = 0) override;

        /// @brief
        /// @return
        [[nodiscard]] size_t multiplicityFactor() const;

        virtual void preProcessWCValues(std::vector<std::complex<double>>& data, bool reverse = false) const;

    protected:

        /// @brief defines new settings for this class
        void defineSettings() override = 0;

        /// @brief adds the index labels for the amplitude tensor for a specific signature to the
        ///        index labels signature list. It will be selected by calling setSignatureIndex
        /// @param[in] tensor  the tensor indices labels
        void addTensor(Tensor&& tensor);

        virtual void updateWilsonCoeffLabelPrefix();

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        std::vector<std::string> _mWCNames;
        std::string _mWCPrefix;
        IndexLabel _mWCLabel{NONE};
        std::vector<Tensor> _tensorList; ///< list of (list of) labels for the tensor indices (one for each signature)
        size_t _multiplicity{1ul};
    };

} // namespace Hammer

#endif
