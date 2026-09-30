///
/// @file  Histos.hh
/// @brief Hammer histogram manager
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_HISTOS_HH
#define HAMMER_HISTOS_HH

#include <cstdint>
#include <vector>

#include "Hammer/Histos.fhh"
#include "Hammer/Config/HammerConfig.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/HammerRoot.fhh"
#include "Hammer/Tools/IOTypes.hh"

#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Math/HistogramDefinition.hh"
#include "Hammer/Math/HistogramSet.hh"

namespace Hammer {

    class Log;
    class DictionaryManager;

    template <typename T>
    using HistoNameDict = std::map<std::string, T>;

    namespace MD = MultiDimensional;

    /// @brief Hammer histogram manager class
    ///
    /// Stores Hammer histograms, ...
    ///
    /// @ingroup Core
    class Histos : public SettingsConsumer {

    public:

        explicit Histos(DictionaryManager* dict = nullptr);

        Histos(const Histos& other) = delete;
        Histos& operator=(const Histos& other) = delete;
        Histos(Histos&& other) = delete;
        Histos& operator=(Histos&& other) = delete;

        ~Histos() noexcept override;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] binSizes
        /// @param[in] hasUnderOverFlow
        /// @param[in] ranges
        /// @param[in] shouldCompress
        /// @param[in] withErrors
        void addHistogramDefinition(const std::string& histogramName, const IndexList& binSizes,
                                    bool hasUnderOverFlow = true, const MD::BinRangeList& ranges = {},
                                    bool shouldCompress = false, bool withErrors = false);
        void addHistogramDefinition(const std::string& histogramName, const MD::BinEdgeList& binEdges,
                                    bool hasUnderOverFlow = true, bool shouldCompress = false, bool withErrors = false);

        void createProjectedHistogram(const std::string& oldName, const std::string& newName,
                                      const std::set<uint16_t>& collapsedIndexPositions);

        void setHistogramCompression(const std::string& histogramName, bool value = true);

        void setHistogramKeepErrors(const std::string& histogramName, bool value = true);

        void addHistogramFixedData(const std::string& histogramName, const MD::SharedTensorData& data);

        void resetHistogramFixedData(const std::string& histogramName);

        void removeHistogram(const std::string& histogramName);

        /// @brief
        void init();

        /// @brief
        /// @param[in] eventId
        void setEventId(EventUID eventId);

        /// @brief
        void clear();

        /// @brief
        /// @param[in] histogramName
        /// @param[in] value
        /// @return
        [[nodiscard]] IndexList getBinIndices(const std::string& histogramName, const MD::BinValue& value) const;

        [[nodiscard]] EventUIDGroup getHistogramEventIds(const std::string& name, const std::string& specialization,
                                                         const std::string& scheme) const;
        [[nodiscard]] MD::BinEdgeList getHistogramEdges(const std::string& histogramName) const;
        [[nodiscard]] IndexList getHistogramShape(const std::string& histogramName) const;

        [[nodiscard]] bool getUnderOverFlows(const std::string& histogramName) const;

        /// @brief
        /// @return
        [[nodiscard]] size_t size() const;

        /// @brief Checks values can be binned
        /// @param[in] name
        /// @param[in] values
        /// @return
        [[nodiscard]] bool canFill(const std::string& name, const std::vector<double>& values) const;

        /// @brief Checks histogram exists
        /// @param[in] name
        /// @return
        void checkExists(const std::string& name) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] IOHistogram getHistogram(const std::string& histogramName, const std::string& schemeName,
                                               const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] EventIdGroupDict<IOHistogram> getHistograms(const std::string& histogramName,
                                                                  const std::string& schemeName,
                                                                  const std::string& specializationName) const;

        [[nodiscard]] bool isValidHistogram(const std::string& histogramName, size_t dim) const;

        [[nodiscard]] std::vector<std::string> getHistogramNames() const;

        [[nodiscard]] bool haveSpecializations(const std::string& name) const;
        [[nodiscard]] std::set<std::string> getSpecializationsForHisto(const std::string& name) const;

        void addSpecializationForHisto(const std::string& histoName, const std::string& specName);

        void removeSpecializationForHisto(const std::string& histoName, const std::string& specName, bool removeData);

        void removeAllSpecializationsForHisto(const std::string& histoName, bool removeData);

#ifdef HAVE_ROOT

        [[nodiscard]] bool isValidHistogram(const std::string& histogramName, const TH1D& h) const;
        [[nodiscard]] bool isValidHistogram(const std::string& histogramName, const TH2D& h) const;
        [[nodiscard]] bool isValidHistogram(const std::string& histogramName, const TH3D& h) const;


        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] std::unique_ptr<TH1D> getHistogram1D(const std::string& histogramName,
                                                           const std::string& schemeName,
                                                           const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] std::unique_ptr<TH2D> getHistogram2D(const std::string& histogramName,
                                                           const std::string& schemeName,
                                                           const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] std::unique_ptr<TH3D> getHistogram3D(const std::string& histogramName,
                                                           const std::string& schemeName,
                                                           const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] EventIdGroupDict<std::unique_ptr<TH1D>>
        getHistograms1D(const std::string& histogramName, const std::string& schemeName,
                        const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] EventIdGroupDict<std::unique_ptr<TH2D>>
        getHistograms2D(const std::string& histogramName, const std::string& schemeName,
                        const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] EventIdGroupDict<std::unique_ptr<TH3D>>
        getHistograms3D(const std::string& histogramName, const std::string& schemeName,
                        const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] rootHistogram
        /// @return
        void setHistogram1D(const std::string& histogramName, const std::string& schemeName, TH1D& rootHistogram,
                            const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] rootHistogram
        /// @return
        void setHistogram2D(const std::string& histogramName, const std::string& schemeName, TH2D& rootHistogram,
                            const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] rootHistogram
        /// @return
        void setHistogram3D(const std::string& histogramName, const std::string& schemeName, TH3D& rootHistogram,
                            const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] rootHistograms
        /// @return
        void setHistograms1D(const std::string& histogramName, const std::string& schemeName,
                             EventIdGroupDict<std::unique_ptr<TH1D>>& rootHistograms,
                             const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] rootHistograms
        /// @return
        void setHistograms2D(const std::string& histogramName, const std::string& schemeName,
                             EventIdGroupDict<std::unique_ptr<TH2D>>& rootHistograms,
                             const std::string& specializationName) const;

        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] rootHistograms
        /// @return
        void setHistograms3D(const std::string& histogramName, const std::string& schemeName,
                             EventIdGroupDict<std::unique_ptr<TH3D>>& rootHistograms,
                             const std::string& specializationName) const;

#endif
        /// @brief
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] binPosition
        /// @param[in] value
        /// @param[in] extraWeight
        void fillHisto(const std::string& histogramName, const std::string& specializationName,
                       const std::string& schemeName, const IndexList& binPosition, Tensor& value,
                       double extraWeight = 1.0);

        //        std::set<EventId> getProcesses(const std::string& histogramName, const std::string& schemeName) const;

    protected:

        /// @brief    logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        /// @brief
        void defineSettings() override;

    public:

        /// @brief
        /// @param[out] msgwriter
        /// @param[in] histogramName
        bool writeDefinition(flatbuffers::FlatBufferBuilder* msgwriter, const std::string& histogramName) const;

        /// @brief
        /// @param[out] msgwriter
        /// @param[in] histogramName
        /// @param[in] schemeName
        /// @param[in] eventIDs
        bool writeHistogram(flatbuffers::FlatBufferBuilder* msgwriter, const std::string& histogramName,
                            const std::string& specializationName, const std::string& schemeName,
                            const EventUID& eventIDs) const;

        /// @brief
        /// @param[in] msgreader
        /// @param[in] merge
        std::string readDefinition(const Serial::FBHistoDefinition* msgreader, bool merge);

        /// @brief
        /// @param[in] msgreader
        /// @param[in] merge
        HistoInfo readHistogram(const Serial::FBHistogram* msgreader, bool merge);


        [[nodiscard]] std::vector<EventUID> getEventIDRepsForHisto(const std::string& name,
                                                                   const std::string& specialization,
                                                                   const std::string& scheme) const;


    private:

        [[nodiscard]] const HistogramSet* getEntry(const std::string& histogramName,
                                                   const std::string& specializationName,
                                                   const std::string& schemeName) const;
        HistogramSet* getEntry(const std::string& histogramName, const std::string& specializationName,
                               const std::string& schemeName);

        [[nodiscard]] std::vector<Tensor> getExternalData(const std::string& schemeName, std::vector<LabelsList> labels,
                                                          const std::string& specializationName) const;

        EventUID _currentEventId;

        HistoNameDict<HistogramDefinition> _histogramDefs;
        HistoNameDict<SpecializationDict<SchemeDict<HistogramSet>>> _histograms;

        HistoNameDict<std::set<std::string>> _usedSpecializations;
        DictionaryManager* _dictionaries;

        bool _initialized{false};
    };

} // namespace Hammer

#endif
