///
/// @file  FormFactorBase.hh
/// @brief Hammer base form factor class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FORMFACTOR_BASE
#define HAMMER_FORMFACTOR_BASE

#ifndef FORMFACTORBASE_FRIENDS
#define FORMFACTORBASE_FRIENDS typedef bool DummyFriend
#endif

#include <string>
#include <utility>
#include <vector>
#include <memory>

#include "Hammer/Math/Integrator.fhh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Tools/ParticleData.hh"
#include "Hammer/Tools/SettingsConsumer.hh"

namespace Hammer {

    class Log;

    /// @brief Base class for form factors
    ///
    /// Provides the form factor contents as tensor...
    ///
    /// @ingroup Base
    class FormFactorBase : public ParticleData, public SettingsConsumer {

        friend class FFEval;
        FORMFACTORBASE_FRIENDS;

    public:

        FormFactorBase();

        // copy constructible but not copy assignable to support MAKE_CLONE
        FormFactorBase(const FormFactorBase& other);
        FormFactorBase& operator=(const FormFactorBase& other) = delete;

        FormFactorBase(FormFactorBase&& other) = delete;
        FormFactorBase& operator=(FormFactorBase&& other) = delete;
        ~FormFactorBase() noexcept override = default;

        /// @brief method to evaluate the object on a specific particle set
        /// @param[in] parent  the parent Particle
        /// @param[in] daughters the daughters (and grand-daughters, if necessary) Particle list
        /// @param[in] references the parent Particle siblings (necessary e.g. for helicity amplitude phase conventions)
        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) override = 0;

        /// @brief
        /// @param[in] label
        /// @return
        virtual std::unique_ptr<FormFactorBase> clone(const std::string& label) = 0;

        std::vector<double> getErrVectorFromDict(const std::map<std::string, double>& errDict) const;
        std::vector<double> getErrVectorFromSettings(bool useDefault = false) const;

        void updateFFErrSettings(const std::vector<double>& values);
        void updateFFErrSettings(const std::map<std::string, double>& values);

        std::map<std::string, double> retrieveFFErrSettings() const;

        void updateFFErrTensor(const std::vector<double>& values, MultiDimensional::SharedTensorData& data) const;

        /// @brief
        /// @return
        const FFPrefixGroup& getFFErrPrefixGroup() const;

        std::pair<FFPrefixGroup, IndexLabel> getFFErrInfo() const;

        /// @brief initializes the form factor (defines settings associated to this form factor, etc.)
        void init();

        /// @brief returns a reference to itself as a Tensor
        /// @return itself
        Tensor& getTensor();

        /// @brief returns a reference to itself as a Tensor
        /// @return itself
        const Tensor& getTensor() const;

        /// @brief
        /// @param[in] intPoints
        /// @return
        Tensor getFFPSIntegrand(const EvaluationGrid& intPoints);

        /// @brief
        /// @return
        const std::string& group() const;

        /// @brief compute unit rescaling of this class wrt MC units
        void calcUnits();

        void bindErrNames();

    protected:

        /// @brief defines new settings for this class
        void defineSettings() override = 0;

        /// @brief
        /// @return
        void defineAndAddErrSettings(const std::vector<std::string>& names);

        /// @brief
        /// @param[in] point other invariant masses
        /// @param[in] masses
        /// @return
        virtual void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) = 0;

        /// @brief Sets the FF parametrization scheme label. Enforced to be same over the set of signatures.
        /// @param[in] name
        void setGroup(const std::string& name);

        /// @brief Sets the XtoY formal process label. May vary over the set of signatures.
        /// @param[in] name
        void setPrefix(const std::string& name);

        /// @brief set units for the FF class parameters
        /// @param[in] name
        void setUnits(const std::string& name);

        /// @brief adds the index labels for the form factor tensor for a specific signature to the
        ///        index labels signature list. It will be selected by calling setSignatureIndex
        /// @param[in] tensor  the tensor indices labels
        void addTensor(Tensor&& tensor);

        /// @brief derived from ParticleData::addProcessSignature, also populates _mFFInfoList
        /// @param[in] parent  the parent PDG code
        /// @param[in] daughters the list of daughter PDG codes
        void addProcessSignature(PdgId parent, const std::vector<PdgId>& daughters);

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        /// @brief check whether FF parameters have been initialized
        [[nodiscard]] bool isInitialized() const {
            return _initialized;
        }
        /// @brief mark FF parameters as initialized
        void setInitialized() const {
            _initialized = true;
        }
        /// @brief reset initialization flag
        void resetInitialized() const {
            _initialized = false;
        }

        FFPrefixGroup _errPrefixGroup;
        IndexLabel _mFFErrLabel{IndexLabel{}};
        std::vector<std::pair<FFPrefixGroup, IndexLabel>>
            _mFFInfoList; ///< list of {prefix/group struct, uncertainty label} (one for each signature)
        std::vector<std::string>* _mFFErrNames{nullptr};
        std::vector<Tensor> _tensorList; ///< list of (list of) labels for the tensor indices (one for each signature)
        double _units{1.};
        bool _unitsinit = false;

    private:

        mutable bool _initialized = false;
    };


    /// @brief Base class for form factors with only one hadron in initial and final states
    ///
    /// Provides a default specialization of the eval method
    ///
    /// @ingroup Base
    class FF1to1Base : public FormFactorBase {
    public:

        FF1to1Base() = default;

        FF1to1Base(const FF1to1Base& other) = default;
        FF1to1Base& operator=(const FF1to1Base& other) = delete;

        FF1to1Base(FF1to1Base&& other) = delete;
        FF1to1Base& operator=(FF1to1Base&& other) = delete;
        ~FF1to1Base() noexcept override = default;

        /// @brief method to evaluate the object on a specific particle set
        /// @param[in] parent  the parent Particle
        /// @param[in] daughters the daughters (and grand-daughters, if necessary) Particle list
        /// @param[in] references the parent Particle siblings (necessary e.g. for helicity amplitude phase conventions)
        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) override;

        /// @brief
        /// @param[in] label
        /// @return
        std::unique_ptr<FormFactorBase> clone(const std::string& label) override = 0;

    protected:

        /// @brief defines new settings for this class
        void defineSettings() override = 0;

        /// @brief
        /// @param[in] point other invariant masses
        /// @param[in] masses
        /// @return
        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override = 0;

        std::tuple<double, double, double> getParentDaughterHadMasses(const std::vector<double>& masses) const;

        static double getW(double Sqq, double Mparent, double Mdaughter);
    };


} // namespace Hammer

#define MAKE_CLONE(OBJ, LABEL)                                  \
    {                                                           \
        auto tmp = std::make_unique<OBJ>(*this);                \
        tmp->setGroup(tmp->group() + "_" + LABEL);              \
        return std::unique_ptr<FormFactorBase>(std::move(tmp)); \
    }

#endif
