///
/// @file  Event.hh
/// @brief Hammer event class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_EVENT_HH
#define HAMMER_EVENT_HH

#include <string>
#include <vector>

#include "Hammer/Histos.fhh"
#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/ProcManager.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"

namespace Hammer {

    class Log;
    class Process;
    class DictionaryManager;

    /// @brief Event container class
    ///
    /// Contains the decay processes within an event, ...
    ///
    /// @ingroup Core
    class Event : public SettingsConsumer {

    public:

        /// @brief
        /// @param[in] histograms
        /// @param[in] dict
        explicit Event(Histos* histograms = nullptr, DictionaryManager* dict = nullptr);

        Event(const Event& other) = delete;
        Event& operator=(const Event& other) = delete;
        Event(Event&& other) = delete;
        Event& operator=(Event&& other) = delete;

        ~Event() noexcept override;

        /// @brief
        /// @param[in] p
        /// @return
        ProcessUID addProcess(Process& p);

        /// @brief
        /// @param[in] id
        void removeProcess(ProcessUID id);

        void setEventBaseWeight(double weight);

        [[nodiscard]] double getEventBaseWeight() const;

        /// @brief
        void clear();

        /// @brief
        void init();

        [[nodiscard]] const EventUID& getEventId() const;

        /// @brief
        /// @param[in] process
        /// @param[in] what
        /// @return
        [[nodiscard]] const Tensor& getAmplitude(ProcessUID process, WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] what
        /// @return
        [[nodiscard]] ProcIdDict<std::reference_wrapper<const Tensor>>
        getAmplitudes(WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] process
        /// @param[in] what
        /// @return
        [[nodiscard]] const Tensor& getSquaredAmplitude(ProcessUID process, WTerm what = WTerm::NUMERATOR) const;

        [[nodiscard]] ProcIdDict<std::reference_wrapper<const Tensor>>
        getSquaredAmplitudes(WTerm what = WTerm::NUMERATOR) const;

        /// @brief
        /// @param[in] schemeName
        /// @param[in] process
        /// @return
        [[nodiscard]] std::vector<std::reference_wrapper<const Tensor>>
        getProcessFormFactors(const std::string& schemeName, ProcessUID process) const;

        /// @brief
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] ProcIdDict<std::vector<std::reference_wrapper<const Tensor>>>
        getFormFactors(const std::string& schemeName) const;

        /// @brief
        /// @param[in] schemeName
        /// @param[in] process
        /// @return
        [[nodiscard]] const Tensor& getTensorWeight(const std::string& schemeName, ProcessUID process,
                                                    const std::string& specializationId = Spec::none()) const;

        /// @brief
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] ProcIdDict<std::reference_wrapper<const Tensor>>
        getTensorWeights(const std::string& schemeName, const std::string& specializationId = Spec::none()) const;

        /// @brief
        /// @param[in] schemeName
        /// @param[in] process
        /// @return
        [[nodiscard]] double getWeight(const std::string& schemeName, ProcessUID process,
                                       const std::string& specializationId = Spec::none()) const;

        /// @brief
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] ProcIdDict<double> getWeights(const std::string& schemeName,
                                                    const std::string& specializationId = Spec::none()) const;

        [[nodiscard]] bool hasWeights() const;

        void removeAllSpecializedWeights();

        void removeSpecializedWeights(const std::string& specializationId);

        /// @brief
        void calc();

        /// @brief
        /// @param[in] label
        /// @param[in] indices
        void setHistogramBin(const std::string& label, const IndexList& indices);

        /// @brief
        void fillHistograms();

        /// @brief
        /// @param[out] msgwriter
        void write(flatbuffers::FlatBufferBuilder* msgwriter) const;

        /// @brief
        /// @param[in] msgreader
        /// @param[in] merge
        bool read(const Serial::FBEvent* msgreader, bool merge);

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        /// @brief
        void defineSettings() override;

        /// @brief
        /// @param[in] id
        /// @return
        [[nodiscard]] const ProcManager& getProcManager(ProcessUID id) const;

        /// @brief
        /// @param[in] id
        /// @return
        ProcManager& getProcManager(ProcessUID id);

        ProcessUID addProcManager(ProcManager&& pm);

    private:

        ProcIdDict<ProcManager> _processes;
        EventUID _eventId;

        SchemeNameList _schemes;

        HistoNameDict<IndexList> _histogramBins;
        Histos* _histograms;
        DictionaryManager* _dictionaries;

        double _eventWeight{1.0};
    };

} // namespace Hammer

#endif
