///
/// @file  ProcManager.hh
/// @brief Container class for all process related data structures
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROCESSMANAGER_HH
#define HAMMER_PROCESSMANAGER_HH

#include <string>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/ProcGraph.hh"
#include "Hammer/ProcRequirements.hh"
#include "Hammer/ProcResults.hh"
#include "Hammer/Process.hh"

namespace Hammer {

    class Log;
    class DictionaryManager;
    class ProcGraph;
    class ProcRequirements;
    class ProcResults;
    class Process;

    /// @brief Decay process class
    ///
    /// Contains the amplitudes, weights and info associated to a decay, ...
    ///
    /// @ingroup Core
    class ProcManager : public SettingsConsumer {

    public:

        ProcManager();

        explicit ProcManager(const Process& inputs);

        explicit ProcManager(const Serial::FBProcess* msgreader);

        ProcManager(const ProcManager& other) = delete;
        ProcManager& operator=(const ProcManager& other) = delete;
        ProcManager(ProcManager&& other) noexcept;
        ProcManager& operator=(ProcManager&& other) = delete;

        ~ProcManager() noexcept override = default;

        /// @brief
        /// @param[in] dictionaries
        virtual void initialize(DictionaryManager* dictionaries);

        void setSettingsHandler(SettingsHandler& sh) override;
        void setSettingsHandler(const SettingsConsumer& sh) override;

        /// @brief
        void calc();

        [[nodiscard]] const ProcResults& results() const;

        ProcResults& modifiableRsults();

        [[nodiscard]] const Process& inputs() const;

        // NumDenPair<double> getPSRates() const;

    protected:

        /// @brief
        /// @param[in] parent
        void calcRates(const ParticleIndex& parent) const;

        /// @brief
        /// @param[in] parent
        [[nodiscard]] Tensor calcPSTensor(const ParticleIndex& parent) const;

        /// @brief
        /// @param[in] parent
        /// @param[in] ffmaps
        /// @param[in] schemeName
        [[nodiscard]] Tensor calcRateTensor(const ParticleIndex& parent, const ProcIdDict<FFIndex>& ffmaps,
                                            const std::string& schemeName) const;
        /// @brief
        /// @param[in] specialization
        /// @param[in] schemeName
        /// @param[in] rate tensor
        [[nodiscard]] Tensor calcSpecRateTensor(const std::string& spec, const std::string& scheme,
                                                const Tensor& rateten) const;

    public:

        /// @brief
        /// @param[in] msgwriter
        /// @param[in] msg
        void write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBProcess>* msg) const;

        /// @brief
        /// @param[in] msgreader
        /// @param[in] merge
        bool read(const Serial::FBProcess* msgreader, bool merge);

    protected:

        /// @brief purely virtual function for a class to define new settings
        void defineSettings() override;

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        void setDictionary(DictionaryManager* dict);

    private:

        std::unique_ptr<Process> _inputs;
        std::unique_ptr<ProcGraph> _graph;
        std::unique_ptr<ProcRequirements> _reqs;
        std::unique_ptr<ProcResults> _results;
        SpecializationDict<SchemeDict<Tensor>>* _processRates;
        //        NumDenPair<double> _processPhaseSpaceRates;

        DictionaryManager* _dictionaries;
    };

} // namespace Hammer

#endif
