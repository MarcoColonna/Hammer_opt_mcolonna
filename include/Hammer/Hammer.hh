///
/// @file  Hammer.hh
/// @brief Main Hammer class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_HAMMER
#define HAMMER_HAMMER

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Hammer/Config/HammerConfig.hh"
#include "Hammer/Tools/HammerRoot.fhh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/Pdg.fhh"
#include "Hammer/Tools/IOTypes.hh"

#include "Hammer/Math/MultiDimensional.fhh"

namespace Hammer {

    class Log;
    class DictionaryManager;
    class Event;
    class SettingsHandler;
    class Histos;
    class Process;
    class FFEval;

    enum class PAction : std::uint8_t { ALL, WEIGHTS, HISTOGRAMS };

    /// @brief Main class
    ///
    /// Contains ...
    ///
    /// @ingroup Core
    class Hammer final : private SettingsConsumer {

    public:

        Hammer();

        Hammer(const Hammer& other) = delete;
        Hammer& operator=(const Hammer& other) = delete;
        Hammer(Hammer&& other) = delete;
        Hammer& operator=(Hammer&& other) = delete;

        ~Hammer() noexcept final;


        /// @brief
        void initRun();

        /// @brief Clears the _event container
        /// @param[in] weight
        void initEvent(double weight = 1.0);

        /// @brief Adds a process to the _event container
        /// @param[in] p
        /// @return HashId of the process
        size_t addProcess(Process& p);


        /// @brief Removes a process to the _event container
        /// @param[in] id
        void removeProcess(size_t id);

        /// @brief
        /// @param[in] name
        /// @param[in] bins
        void setEventHistogramBin(const std::string& name, const IndexList& bins);

        void fillEventHistogram(const std::string& name, const std::vector<double>& values);

        /// @brief
        /// @param[in] weight
        void setEventBaseWeight(double weight);

        /// @brief
        void processEvent(PAction what = PAction::ALL);


        /// @brief
        /// @param[in] buffer
        /// @param[in] merge
        bool loadEventWeights(const IOBuffer& buffer, bool merge = false);

        /// @brief
        /// @return
        [[nodiscard]] IOBuffer saveEventWeights() const;

        /// @brief
        /// @param[in] buffer
        /// @param[in] merge
        bool loadRunHeader(IOBuffer& buffer, bool merge = false);

        /// @brief
        /// @return
        [[nodiscard]] IOBuffer saveRunHeader() const;

        /// @brief
        /// @param[in] buffer
        /// @param[in] merge
        /// @return
        std::string loadHistogramDefinition(IOBuffer& buffer, bool merge = false);

        /// @brief
        /// @param[in] buffer
        /// @param[in] merge
        HistoInfo loadHistogram(IOBuffer& buffer, bool merge = false);

        /// @brief
        /// @param[in] name
        /// @return
        [[nodiscard]] IOBuffers saveHistogram(const std::string& name) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @return
        [[nodiscard]] IOBuffers saveHistogram(const std::string& name, const std::string& scheme) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @return
        [[nodiscard]] IOBuffers saveHistogram(const std::string& name, const std::string& scheme,
                                              const std::string& specialization) const;

        /// @brief
        /// @param[in] name
        /// @param[in] eventIDs
        /// @return
        [[nodiscard]] IOBuffers saveHistogram(const std::string& name, const EventUIDGroup& eventIDs) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] eventIDs
        /// @return
        [[nodiscard]] IOBuffers saveHistogram(const std::string& name, const std::string& scheme,
                                              const EventUIDGroup& eventIDs) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] eventIDs
        /// @return
        [[nodiscard]] IOBuffers saveHistogram(const std::string& name, const std::string& scheme,
                                              const std::string& specialization, const EventUIDGroup& eventIDs) const;

        /// @brief
        /// @param[in] info
        /// @return
        [[nodiscard]] IOBuffers saveHistogram(const HistoInfo& info) const;

        /// @brief
        /// @param[in] buffer
        /// @param[in] merge
        bool loadRates(IOBuffer& buffer, bool merge = false);

        /// @brief
        /// @return
        [[nodiscard]] IOBuffer saveRates() const;

        /// @brief Reads settings in the cards
        /// @param[in] fileDecays
        /// @param[in] fileOptions
        void readCards(const std::string& fileDecays, const std::string& fileOptions);

        /// @brief Reads settings in the cards
        /// @param[in] fileOptions
        /// @param[in] useDefault
        void saveOptionCard(const std::string& fileOptions, bool useDefault = false) const;

        /// @brief Reads settings in the cards
        /// @param[in] fileDecays
        void saveHeaderCard(const std::string& fileDecays) const;

        void saveReferences(const std::string& fileRefs) const;

        /// @brief
        /// @param[in] options
        void setOptions(const std::string& options);

        /// @brief
        /// @param[in] options
        void setHeader(const std::string& options);

        /// @brief Add total weight sum histogram with compression and errors
        /// @details Histogram is uniquely identified by "Total Sum of Weights", and is a (single bin) zero dim
        /// histogram
        /// @param[in] compress
        /// @param[in] witherrors
        void addTotalSumOfWeights(bool compress = false, bool witherrors = false);

        /// @brief Adds a tensor histogram
        /// @details Histogram is uniquely identified by a name. The dimensionality is arbitrary, with bin entries
        /// generically of type Tensor Example: hammer.addHistogram("q2VsEmuon", {10,20}, false, {{0.,12.},{0.,2.4}})
        /// @param[in] name
        /// @param[in] binSizes
        /// @param[in] hasUnderOverFlow
        /// @param[in] ranges
        void addHistogram(const std::string& name, const IndexList& binSizes, bool hasUnderOverFlow = true,
                          const std::vector<std::pair<double, double>>& ranges = {});

        /// @brief Adds a tensor histogram
        /// @details Histogram is uniquely identified by a name. The dimensionality is arbitrary, with bin entries
        /// generically of type Tensor Example: hammer.addHistogram("q2VsEmuon", {10,20}, false)
        /// @param[in] name
        /// @param[in] binEdges
        /// @param[in] hasUnderOverFlow
        void addHistogram(const std::string& name, const std::vector<std::vector<double>>& binEdges,
                          bool hasUnderOverFlow);

        void collapseProcessesInHistogram(const std::string& name);

        void keepErrorsInHistogram(const std::string& name, bool value = true);

        void reconcileSpecializations() const;

        void saveEventGeneralWeights(bool save = true);

        void reconcileGeneralWeightSetting() const;

        void reconcileGeneralHistogramSetting(const std::string& histogramName) const;

        void createWCSpecialization(const std::string& name, const std::string& wcspace,
                                    const std::vector<std::string>& coordinates = {}) const;

        void setWCSpecializationOrigin(const std::string& name, const std::string& wcspace,
                                       const std::map<std::string, std::complex<double>>& settings) const;

        void setWCSpecializationOrigin(const std::string& name, const std::string& wcspace,
                                       const std::vector<std::complex<double>>& values) const;

        void setWCSpecializationCoord(const std::string& name, const std::string& wcspace, const std::string& coord,
                                      const std::map<std::string, std::complex<double>>& settings) const;

        void setWCSpecializationCoord(const std::string& name, const std::string& wcspace, const std::string& coord,
                                      const std::vector<std::complex<double>>& values) const;

        void setWCSpecializationBasis(const std::string& name, const std::string& wcspace,
                                      const std::vector<std::map<std::string, std::complex<double>>>& subspace) const;

        void removeWCSpecialization(const std::string& name, bool removeExistingData = false) const;

        void removeAllWCSpecializations(bool removeExistingData = false) const;

        [[nodiscard]] std::set<std::string> availableWCSpecializationIds() const;


        void applyWCSpecializationInWeights(const std::string& name) const;

        void removeWCSpecializationInWeights(const std::string& name) const;

        [[nodiscard]] std::set<std::string> appliedWCSpecializationsInWeights() const;


        void applyWCSpecializationInHistogram(const std::string& histogramName, const std::string& specName) const;

        void applyWCSpecializationInAllHistograms(const std::string& specName) const;

        void removeWCSpecializationInHistogram(const std::string& histogramName, const std::string& specName) const;

        [[nodiscard]] std::set<std::string> appliedWCSpecializationsInHistogram(const std::string& histogramName) const;


        void specializeFFInHistogram(const std::string& name, const std::string& process, const std::string& group,
                                     const std::vector<double>& values);

        void specializeFFInHistogram(const std::string& name, const std::string& process, const std::string& group,
                                     const std::map<std::string, double>& settings);


        void removeFFSpecializationInHistogram(const std::string& name);

        void createProjectedHistogram(const std::string& oldName, const std::string& newName,
                                      const std::set<uint16_t>& collapsedIndexPositions);

        void removeHistogram(const std::string& name);

        /// @brief Shows all known FF classes
        /// @details Lists by group name for all class prefixes
        void showAvailableFFParams() const;

        void showAvailableFFParams(const std::string& prefix) const;

        /// @brief Adds a form factor scheme
        /// @details A form factor scheme is a name (string) plus a map from a subprocess name to a FF parametrization
        /// Example: hammer.addFFScheme("MixedBag", {{"BD", "CLN"}, {"BD*", "ISGW2"}})
        /// @param[in] schemeName
        /// @param[in] schemes
        void addFFScheme(const std::string& schemeName, const std::map<std::string, std::string>& schemes);

        /// @brief Sets the FF schemes for the denominator
        /// @details A map from each hadronic subprocess to a specified FF parametrization
        /// Example: hammer.setFFInputScheme({{"B0barD+", "ISGW2"}})
        /// @param[in] schemes
        void setFFInputScheme(const std::map<std::string, std::string>& schemes);

        /// @brief Removes a form factor scheme
        /// Example: hammer.removeFFScheme("MixedBag")
        /// @param[in] schemeName
        void removeFFScheme(const std::string& schemeName);

        [[nodiscard]] std::vector<std::string> getFFSchemeNames() const;

        /// @brief Adds a combinatoric set to included process specifications
        /// @param[in] names
        void includeDecay(const std::vector<std::string>& names);

        /// @brief Adds an included process specification
        /// @param[in] name
        void includeDecay(const std::string& name);

        /// @brief Adds a combinatoric set to forbidden process specifications
        /// @param[in] names
        void forbidDecay(const std::vector<std::string>& names);

        /// @brief Adds a forbidden process specification
        /// @param[in] name
        void forbidDecay(const std::string& name);

        /// @brief Sets pure PS by vertex
        /// @details A map from each vertex to specified position: "NUMERATOR", "DENOMINATOR", "COMMON"
        /// @param[in] vertices
        /// @param[in] what
        void addPurePSVertices(const std::set<std::string>& vertices, WTerm what = WTerm::NUMERATOR);

        void clearPurePSVertices(WTerm what = WTerm::NUMERATOR);

        /// @brief Sets the units
        /// @details default is GeV
        /// @param[in] name
        void setUnits(const std::string& name = "GeV");

        /// @brief
        /// @param[in] process
        /// @param[in] group
        /// @param[in] names
        void renameFFEigenvectors(const std::string& process, const std::string& group,
                                  const std::vector<std::string>& names);

        /// @brief
        /// @param[in] wcspace
        /// @param[in] values
        /// @param[in] what
        void setWilsonCoefficients(const std::string& wcspace, const std::vector<std::complex<double>>& values,
                                   WTerm what = WTerm::NUMERATOR);

        /// @brief
        /// @param[in] wcspace
        /// @param[in] settings
        /// @param[in] what
        void setWilsonCoefficients(const std::string& wcspace,
                                   const std::map<std::string, std::complex<double>>& settings,
                                   WTerm what = WTerm::NUMERATOR);

        /// @brief
        /// @param[in] wcspace
        /// @param[in] values
        void setWilsonCoefficientsLocal(const std::string& wcspace, const std::vector<std::complex<double>>& values);

        /// @brief
        /// @param[in] wcspace
        /// @param[in] settings
        void setWilsonCoefficientsLocal(const std::string& wcspace,
                                        const std::map<std::string, std::complex<double>>& settings);

        /// @brief
        /// @param[in] wcspace
        /// @param[in] what
        void resetWilsonCoefficients(const std::string& wcspace, WTerm what = WTerm::NUMERATOR);

        [[nodiscard]] std::map<std::string, std::complex<double>>
        retrieveWilsonCoefficients(const std::string& wcspace, WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] process
        /// @param[in] group
        /// @param[in] values
        void setFFEigenvectors(const std::string& process, const std::string& group, const std::vector<double>& values);

        /// @brief
        /// @param[in] process
        /// @param[in] group
        /// @param[in] settings
        void setFFEigenvectors(const std::string& process, const std::string& group,
                               const std::map<std::string, double>& settings);

        /// @brief
        /// @param[in] process
        /// @param[in] group
        /// @param[in] values
        void setFFEigenvectorsLocal(const std::string& process, const std::string& group,
                                    const std::vector<double>& values);

        /// @brief
        /// @param[in] process
        /// @param[in] group
        /// @param[in] settings
        void setFFEigenvectorsLocal(const std::string& process, const std::string& group,
                                    const std::map<std::string, double>& settings);

        /// @brief
        /// @param[in] process
        /// @param[in] group
        void resetFFEigenvectors(const std::string& process, const std::string& group);

        [[nodiscard]] std::map<std::string, double> retrieveFFEigenvectors(const std::string& process,
                                                                           const std::string& group) const;

        /// @brief
        /// @param[in] scheme
        /// @param[in] processes
        /// @param[in] specialization
        /// @return
        [[nodiscard]] double getWeight(const std::string& scheme, const std::vector<HashId>& processes = {},
                                       const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] scheme
        /// @param[in] processes
        /// @param[in] specialization
        /// @return
        [[nodiscard]] double getWeight(const std::string& scheme, const std::string& specialization) const;

        /// @brief
        /// @param[in] scheme
        /// @param[in] processes
        /// @param[in] specialization
        /// @return
        [[nodiscard]] double getWeight(const std::string& scheme,
                                       const std::vector<std::vector<std::string>>& processes,
                                       const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] std::map<size_t, double> getWeights(const std::string& scheme,
                                                          const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] id
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return rate double
        [[nodiscard]] double getRate(const HashId& id, const std::string& scheme,
                                     const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] parent pdg id
        /// @param[in] daughters pdg ids
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return rate double
        [[nodiscard]] double getRate(const PdgId& parent, const std::vector<PdgId>& daughters,
                                     const std::string& scheme, const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] vertex string
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return rate double
        [[nodiscard]] double getRate(const std::string& vertex, const std::string& scheme,
                                     const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] id
        /// @param[in] specialization
        /// @return rate double
        [[nodiscard]] double getDenominatorRate(const HashId& id) const;

        /// @brief
        /// @param[in] parent pdg id
        /// @param[in] daughters pdg ids
        /// @return rate double
        [[nodiscard]] double getDenominatorRate(const PdgId& parent, const std::vector<PdgId>& daughters) const;

        /// @brief
        /// @param[in] vertex string
        /// @return rate double
        [[nodiscard]] double getDenominatorRate(const std::string& vertex) const;

        /// @brief
        /// @param[in] scheme
        /// @param[in] process
        /// @param[in] PS point
        /// @param[in] masses
        /// @return vector of FF doubles
        [[nodiscard]] std::vector<double> evalFormFactors(const std::string& scheme, const std::string& process, 
                                                           const std::vector<double>& point, const std::vector<double>& masses = {}) const;
        
        /// @brief
        /// @param[in] scheme
        /// @param[in] parent
        /// @param[in] daughters
        /// @param[in] PS point
        /// @param[in] masses
        /// @return vector of FF doubles
        [[nodiscard]] std::vector<double> evalFormFactors(const std::string& scheme, const PdgId& parent, 
                                                                const std::vector<PdgId>& daughters, const std::vector<double>& point, 
                                                                const std::vector<double>& masses = {}) const;
        
        
        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] IOHistogram getHistogram(const std::string& name, const std::string& scheme,
                                               const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] EventIdGroupDict<IOHistogram>
        getHistograms(const std::string& name, const std::string& scheme,
                      const std::string& specialization = Spec::none()) const;

        [[nodiscard]] EventUIDGroup getHistogramEventIds(const std::string& name, const std::string& scheme,
                                                         const std::string& specialization = Spec::none()) const;
        [[nodiscard]] std::vector<std::vector<double>> getHistogramBinEdges(const std::string& name) const;
        [[nodiscard]] IndexList getHistogramShape(const std::string& name) const;
        [[nodiscard]] bool histogramHasUnderOverFlows(const std::string& name) const;

#ifdef HAVE_ROOT

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] std::unique_ptr<TH1D> getHistogram1D(const std::string& name, const std::string& scheme,
                                                           const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] std::unique_ptr<TH2D> getHistogram2D(const std::string& name, const std::string& scheme,
                                                           const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] std::unique_ptr<TH3D> getHistogram3D(const std::string& name, const std::string& scheme,
                                                           const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] EventIdGroupDict<std::unique_ptr<TH1D>>
        getHistograms1D(const std::string& name, const std::string& scheme,
                        const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] EventIdGroupDict<std::unique_ptr<TH2D>>
        getHistograms2D(const std::string& name, const std::string& scheme,
                        const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @return
        [[nodiscard]] EventIdGroupDict<std::unique_ptr<TH3D>>
        getHistograms3D(const std::string& name, const std::string& scheme,
                        const std::string& specialization = Spec::none()) const;


        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @param[out] histogram
        /// @return
        void setHistogram1D(const std::string& name, const std::string& scheme, TH1D& histogram,
                            const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @param[out] histogram
        /// @return
        void setHistogram2D(const std::string& name, const std::string& scheme, TH2D& histogram,
                            const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @param[out] histogram
        /// @return
        void setHistogram3D(const std::string& name, const std::string& scheme, TH3D& histogram,
                            const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @param[out] histograms
        /// @return
        void setHistograms1D(const std::string& name, const std::string& scheme,
                             EventIdGroupDict<std::unique_ptr<TH1D>>& histograms,
                             const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @param[out] histograms
        /// @return
        void setHistograms2D(const std::string& name, const std::string& scheme,
                             EventIdGroupDict<std::unique_ptr<TH2D>>& histograms,
                             const std::string& specialization = Spec::none()) const;

        /// @brief
        /// @param[in] name
        /// @param[in] scheme
        /// @param[in] specialization
        /// @param[out] histograms
        /// @return
        void setHistograms3D(const std::string& name, const std::string& scheme,
                             EventIdGroupDict<std::unique_ptr<TH3D>>& histograms,
                             const std::string& specialization = Spec::none()) const;


#endif

    private:

        /// @brief purely virtual function for a class to define new settings
        void defineSettings() final;

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        std::unique_ptr<DictionaryManager> _containers;
        std::unique_ptr<SettingsHandler> _settings;
        std::unique_ptr<Histos> _histograms;
        std::unique_ptr<flatbuffers::FlatBufferBuilder> _builder;
        std::unique_ptr<Event> _event;
        std::unique_ptr<FFEval> _ffeval;
        double _mcunits = 1.;
        bool _postInitRun{false};
        bool _postLoadRunHeader{false};
    };

} // namespace Hammer

#endif
