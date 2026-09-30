///
/// @file  RateBase.hh
/// @brief Hammer base rate class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_BASE
#define HAMMER_RATE_BASE

#ifndef RATEBASE_FRIENDS
#define RATEBASE_FRIENDS typedef bool DummyFriend
#endif

#include <string>
#include <vector>

#include "Hammer/Particle.fhh"
#include "Hammer/Math/Integrator.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Tools/ParticleData.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/Logging.hh"

namespace Hammer {

    class Integrator;

    /// @brief Base class for rates
    ///
    /// Provides the rate of a decay process as a N+1-dimensional tensor, whose last index is the
    /// is the \f$ d\Gamma/dq^2 \f$ (up to wilson coefficients and form factors) evaluated at different \f$ q^2 \f$
    /// values, so that the \f$ q^2 \f$ integration can be expressed as \f$ \mbox{Tr} (R \cdot w \cdot FF) \f$,
    /// where \f$ R \f$ is the rate tensor, \f$ FF \f$ the (outer) squared form factor tensor evaluated at the same \f$
    /// q^2 \f$ points and \f$ W \f$ is a diagonal matrix encoding the weights for the gaussian quadrature integration.
    /// Internally, rate units are always in GeV (same as Pdg class partial widths etc).
    ///
    /// @ingroup Base
    class RateBase : public ParticleData, public SettingsConsumer {

        RATEBASE_FRIENDS;

    public:

        RateBase();

        RateBase(const RateBase& other) = delete;
        RateBase& operator=(const RateBase& other) = delete;

        RateBase(RateBase&& other) = delete;
        RateBase& operator=(RateBase&& other) = delete;
        ~RateBase() noexcept override = default;

    private:

        /// @brief method to evaluate the object on a specific particle set. Implementation necessary from ParticleData,
        ///        but empty for rates because they are pre-computed in the init method
        /// @param[in] parent  the parent Particle
        /// @param[in] daughters the daughters (and grand-daughters, if necessary) Particle list
        /// @param[in] references the parent Particle siblings (necessary e.g. for helicity amplitude phase conventions)
        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) override;

    public:

        /// @brief initializes the rate (defines settings associated to this rate, computes the rate tensor)
        void init();

        /// @brief returns string of pdg ids in the process
        /// @return string
        [[nodiscard]] std::string getPdgIdString() const;

        /// @brief returns pdg ids in the process
        /// @return parent Id, daughter Ids
        [[nodiscard]] std::pair<PdgId, std::vector<PdgId>> getPdgIds() const;

        /// @brief returns a reference to itself as a Tensor
        /// @return itself
        Tensor& getTensor();

        /// @brief returns a reference to itself as a Tensor
        /// @return itself
        [[nodiscard]] const Tensor& getTensor() const;

        /// @brief return the physical range of \f$ Q^2 \f$ for this process (to be used set the integration limits)
        /// @return \f$ [Q^2_{min}, Q^2_{max}] \f$
        [[nodiscard]] const IntegrationBoundaries& getIntegrationBoundaries() const;

        [[nodiscard]] const EvaluationGrid& getEvaluationPoints() const;

        /// @brief evaluates the rate by computing the (rank N) rate tensor at different \f$ q^2 \f$ points
        ///        using evalAtPSPoint on the integration points determined by Integrator and creating a rank N+1
        ///        tensor.
        void calcTensor();

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        /// @brief defines new settings for this class
        void defineSettings() override;

        /// @brief evaluates the rate at a specific point in \f$ q^2 \f$ as a rank N Tensor
        /// @param[in] point  the \f$ q^2 \f$ value
        /// @return the rate tensor
        virtual Tensor evalAtPSPoint(const std::vector<double>& point) = 0;

        /// @brief adds the \f$ q^2 \f$ integration limits for a specific signature to the
        ///        \f$ q^2 \f$ range list. It will be selected by calling setSignatureIndex
        /// @param[in] boundaries  \f$ q_{min} \f$
        void addIntegrationBoundaries(const IntegrationBoundaries& boundaries);

        /// @brief adds the index labels for the amplitude tensor for a specific signature to the
        ///        index labels signature list. It will be selected by calling setSignatureIndex
        /// @param[in] tensor  the tensor indices labels
        void addTensor(Tensor&& tensor);

        Integrator _integ;
        std::vector<IntegrationBoundaries>
            _mPSRanges;                  ///< list of \f$ q^2 \f$ integration limits (one for each signature)
        IndexType _nPoints;              ///< number of \f$ q^2 \f$ points used for the numerical integration
        std::vector<Tensor> _tensorList; ///< list of (list of) labels for the tensor indices (one for each signature)
    };

    namespace PS {

        inline BoundaryFunction makeQ2Function(double qmin, double qmax) {
            return ([qmin, qmax](const std::vector<double>&) -> std::pair<double, double> {
                return std::make_pair(pow(qmin, 2.), pow(qmax, 2.));
            });
        }
    } // namespace PS

} // namespace Hammer

#endif
