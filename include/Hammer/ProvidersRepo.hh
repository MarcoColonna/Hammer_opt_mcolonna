///
/// @file  ProvidersRepo.hh
/// @brief Interface class for amplitudes, rates, FFs dictionary container
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROVIDERSREPO_HH
#define HAMMER_PROVIDERSREPO_HH

#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/Pdg.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"

namespace Hammer {

    class AmplitudeBase;
    class FormFactorBase;
    class RateBase;
    class SchemeDefinitions;

    class IWCFFErrProviders {
    public:

        virtual ~IWCFFErrProviders() = default;

        /// @brief
        /// @param[in] parent
        /// @param[in] daughters
        /// @param[in] granddaughters
        /// @return
        [[nodiscard]] virtual AmplitudeBase* getAmplitude(PdgId parent, const std::vector<PdgId>& daughters,
                                                          const std::vector<PdgId>& granddaughters = {}) const = 0;

        [[nodiscard]] virtual AmplitudeBase* getWCProvider(const std::string& wcPrefix) const = 0;
        [[nodiscard]] virtual IndexLabel getWCLabel(const std::string& wcPrefix) const = 0;
        [[nodiscard]] virtual std::map<IndexLabel, AmplitudeBase*> getAllWCProviders() const = 0;

        [[nodiscard]] virtual FormFactorBase* getFFErrProvider(const FFPrefixGroup& process) const = 0;
        [[nodiscard]] virtual IndexLabel getFFErrLabel(const FFPrefixGroup& process) const = 0;
        [[nodiscard]] virtual std::map<IndexLabel, SchemeDict<FormFactorBase*>> getAllFFErrProviders() const = 0;

        [[nodiscard]] virtual std::set<std::string> schemeNamesFromPrefixAndGroup(const FFPrefixGroup& value) const = 0;
    };

    /// @brief Main class
    ///
    /// Contains ...
    ///
    /// @ingroup Core
    class ProvidersRepo : public IWCFFErrProviders, public SettingsConsumer {

    public:

        explicit ProvidersRepo(const SchemeDefinitions* schemeDefs);

        ProvidersRepo(const ProvidersRepo& other) = delete;
        ProvidersRepo& operator=(const ProvidersRepo& other) = delete;
        ProvidersRepo(ProvidersRepo&& other) = delete;
        ProvidersRepo& operator=(ProvidersRepo&& other) = delete;

        ~ProvidersRepo() noexcept override;

        void setSettingsHandler(SettingsHandler& sh) override;

        void initialize();

        /// @brief
        /// @param[in] parent
        /// @param[in] daughters
        /// @param[in] granddaughters
        /// @return
        [[nodiscard]] AmplitudeBase* getAmplitude(PdgId parent, const std::vector<PdgId>& daughters,
                                                  const std::vector<PdgId>& granddaughters = {}) const override;


        /// @brief Obtains a rate class
        /// @details Extracts the pointer to the rate class from the container _rates, corresponding to the HashId of
        /// the parent plus daughters (and granddaughters if any)
        /// @param[in] parent
        /// @param[in] daughters
        /// @param[in] granddaughters
        /// @return Pointer to the rate class
        [[nodiscard]] virtual RateBase* getRate(PdgId parent, const std::vector<PdgId>& daughters,
                                                const std::vector<PdgId>& granddaughters = {}) const;

        /// @brief Obtains a partial width
        /// @details Extracts the pointer to the tensor from the container _partialWidths, corresponding to the HashId
        /// of the parent plus daughters (and granddaughters if any)
        /// @param[in] parent
        /// @param[in] daughters
        /// @param[in] granddaughters
        /// @return pointer to partial width double
        [[nodiscard]] virtual const double* getPartialWidth(PdgId parent, const std::vector<PdgId>& daughters,
                                                            const std::vector<PdgId>& granddaughters = {}) const;

        // /// @brief Obtains FF classes
        // /// @details Extracts the vector of pointers to the chosen FF parametrization classes from the container
        // _formFactors, corresponding to the HashId of the parent plus daughters (and granddaughters)
        // /// @param[in] parent
        // /// @param[in] daughters
        // /// @param[in] granddaughters
        // /// @return Vector of pointers to chosen FF parameterizations
        // std::vector<FormFactorBase*> getFormFactor(PdgId parent, const std::vector<PdgId>& daughters,
        //                                            const std::vector<PdgId>& granddaughters = {}) const;

        /// @brief For each pair of charge conjugate processes, form factor signature dicts contain only one as a key
        /// @param[in] processId
        /// @return 
        [[nodiscard]] bool isHadronicIdKey(const HashId& processId) const;

        /// @brief
        /// @param[in] processId
        /// @return
        [[nodiscard]] virtual std::vector<FormFactorBase*> getFormFactor(HadronicUID processId) const;


        [[nodiscard]] AmplitudeBase* getWCProvider(const std::string& wcPrefix) const override;
        [[nodiscard]] IndexLabel getWCLabel(const std::string& wcPrefix) const override;
        [[nodiscard]] std::map<IndexLabel, AmplitudeBase*> getAllWCProviders() const override;

        [[nodiscard]] FormFactorBase* getFFErrProvider(const FFPrefixGroup& process) const override;
        [[nodiscard]] IndexLabel getFFErrLabel(const FFPrefixGroup& process) const override;
        [[nodiscard]] std::map<IndexLabel, SchemeDict<FormFactorBase*>> getAllFFErrProviders() const override;

        /// @brief Clones a FF parametrization
        /// @param[in] processId
        /// @param[in] baseIndex
        /// @param[in] newLabel
        /// @return
        FFIndex duplicateFormFactor(HadronicUID processId, FFIndex baseIndex, const std::string& newLabel);

        void initPreSchemeDefs();
        void initPostSchemeDefs();

        [[nodiscard]] bool checkFFPrefixAndGroup(const FFPrefixGroup& value) const;

        [[nodiscard]] std::set<std::string> schemeNamesFromPrefixAndGroup(const FFPrefixGroup& value) const override;

        [[nodiscard]] HadronicIdDict<std::vector<std::string>> getFFGroups() const;

        [[nodiscard]] std::map<std::string, std::set<std::string>> getAvailableFFParams() const;

        void showAvailableFFParams() const;

        void showAvailableFFParams(const std::string& prefix) const;

        [[nodiscard]] bool renameFFEigenvectors(const FFPrefixGroup& process,
                                                const std::vector<std::string>& names) const;

    private:

        void fixWCValues();

        void defineSettings() override;

        static Log& getLog();

        void initializeWCs();
        void initializeFFErrs();

        /// @brief Known amplitudes
        /// @details A map {HashId (parent, daughters (and granddaughters)), amplitude class pointer} for all amplitudes
        /// known to hammer
        AmplitudeIdDict<std::unique_ptr<AmplitudeBase>> _amplitudes;
        /// @brief Known rates
        /// @details A map {HashId (parent, daughters), rate class pointer} for all rates known to hammer
        VertexIdDict<std::unique_ptr<RateBase>> _rates;
        /// @brief Known formfactors
        /// @details A map {HashId hadronicId(), vector of FF class pointers} for all FFs known to hammer
        HadronicIdDict<std::vector<std::unique_ptr<FormFactorBase>>> _formFactors;
        /// @brief Known partial widths
        /// @details A map {HashId (parent, daughters), tensor pointer} for all partial widths known to hammer
        VertexIdDict<double> _partialWidths;

        FFPrefixGroupDict<std::set<std::string>> _schemeNamesFromPrefixGroup;

        WCPrefixDict<std::pair<IndexLabel, AmplitudeBase*>> _mWCProviders;
        FFPrefixGroupDict<std::pair<IndexLabel, FormFactorBase*>> _mFFErrProviders;

        const SchemeDefinitions* _schemeDefs;
    };

} // namespace Hammer

#endif
