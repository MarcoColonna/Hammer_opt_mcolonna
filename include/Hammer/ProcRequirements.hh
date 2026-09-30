///
/// @file  ProcRequirements.hh
/// @brief Container class for required ingredients for the process weight calculation
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROCREQUIREMENTS_HH
#define HAMMER_PROCREQUIREMENTS_HH

#include <string>
#include <vector>

#include "Hammer/IndexTypes.hh"
#include "Hammer/ProcGraph.fhh"
#include "Hammer/Math/Tensor.hh"

namespace Hammer {

    class Log;
    class DictionaryManager;
    class AmplitudeBase;
    class RateBase;
    class FormFactorBase;
    class Process;

    /// @brief Decay process class
    ///
    /// Contains the amplitudes, weights and info associated to a decay, ...
    ///
    /// @ingroup Core
    class ProcRequirements {

    public:

        ProcRequirements() = default;

        ProcRequirements(const ProcRequirements& other) = delete;
        ProcRequirements& operator=(const ProcRequirements& other) = delete;
        ProcRequirements(ProcRequirements&& other) = delete;
        ProcRequirements& operator=(ProcRequirements&& other) = delete;

        ~ProcRequirements() noexcept = default;

        /// @brief
        /// @param[in] dictionaries
        /// @param[in] inputs
        /// @param[in] graph
        size_t initialize(const DictionaryManager* dictionaries, const Process* inputs, const ProcGraph* graph);

        /// @brief
        /// @return
        [[nodiscard]] std::vector<
            std::tuple<ParticleIndex, ParticleIndex, NumDenPair<AmplEntry>, NumDenPair<double>, NumDenPair<double>>>
        generatedAmplsMultsPS() const;

        [[nodiscard]] std::pair<double, double> getPSRates() const;

        [[nodiscard]] double calcCorrectionFactor(WTerm what = WTerm::NUMERATOR) const;

        [[nodiscard]] const VertexDict<SelectedAmplEntry>& amplitudes() const;
        [[nodiscard]] const VertexDict<RateBase*>& rates() const;
        [[nodiscard]] const VertexDict<FFIndexDict<FormFactorBase*>>& formFactors() const;
        [[nodiscard]] const std::vector<Tensor>& denominatorWilsonCoeffs() const;
        [[nodiscard]] const std::vector<Tensor>& denominatorFFEigenVectors() const;
        [[nodiscard]] const std::vector<std::reference_wrapper<const Tensor>>&
        specializedWilsonCoeffs(const std::string& id) const;
        [[nodiscard]] const VertexDict<HashId>& rateIds() const;
        [[nodiscard]] const VertexDict<const double*>& partialWidths() const;

        [[nodiscard]] const std::set<std::string>& availableSpecializations() const;

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        ///@brief
        /// @param[in] descendant
        [[nodiscard]] std::pair<ParticleIndex, bool> getAncestorId(ParticleIndex descendant) const;

        /// @brief
        void getAmplitudes();

        /// @brief
        void getFormFactors();

        /// @brief
        void getRates();

    private:

        VertexDict<VertexUID> _rateIds;

        const DictionaryManager* _dictionaries = nullptr;
        const ProcGraph* _graph = nullptr;
        const Process* _inputs = nullptr;

        VertexDict<SelectedAmplEntry> _requiredAmplitudes;
        VertexDict<RateBase*> _requiredRates;
        VertexDict<FFIndexDict<FormFactorBase*>> _requiredFormFactors;
        VertexDict<const double*> _requiredPWs;
        VertexDict<NumDenPair<double>> _multPSFactors;
        VertexDict<NumDenPair<double>> _massPSFactors;
        std::vector<Tensor> _denominatorWilsonCoeffs;
        std::vector<Tensor> _denominatorFFEigenVectors;
        SpecializationDict<std::vector<std::reference_wrapper<const Tensor>>> _specializedWilsonCoeffs;

        std::set<std::string> _usedSpecializations;
    };

} // namespace Hammer

#endif
