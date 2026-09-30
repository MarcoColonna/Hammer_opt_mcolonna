///
/// @file  ProcRates.hh
/// @brief Container class for process rate tensors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_PROCRATES_HH
#define HAMMER_PROCRATES_HH

#include <string>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/HammerSerial.fhh"

namespace Hammer {

    class Log;
    class ExternalData;

    /// @brief Decay process class
    ///
    /// Contains the amplitudes, weights and info associated to a decay, ...
    ///
    /// @ingroup Core
    class ProcRates : public SettingsConsumer {

    public:

        explicit ProcRates(const ExternalData* ext);

        explicit ProcRates(const Serial::FBProcData* msgreader);

        ProcRates(const ProcRates& other) = delete;
        ProcRates& operator=(const ProcRates& other) = delete;
        ProcRates(ProcRates&& other) = delete;
        ProcRates& operator=(ProcRates&& other) = delete;

        ~ProcRates() noexcept override;

        /// @brief
        /// @param[in] id
        /// @return
        virtual SpecializationDict<SchemeDict<Tensor>>* getProcessRates(ProcessUID id);

        /// @brief
        /// @param[in] id
        /// @param[in] schemeName
        /// @return
        [[nodiscard]] virtual double getVertexRate(const ProcessUID& id, const std::string& schemeName,
                                                   const std::string& specializationId = Spec::none()) const;

        void init();

        /// @brief
        /// @param[in] msgwriter
        void write(flatbuffers::FlatBufferBuilder* msgwriter) const;

        /// @brief
        /// @param[in] msgreader
        /// @param[in] merge
        bool read(const Serial::FBRates* msgreader, bool merge);

    protected:

        /// @brief purely virtual function for a class to define new settings
        void defineSettings() override;

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    private:

        ProcIdDict<SpecializationDict<SchemeDict<Tensor>>> _processRates;
        const ExternalData* _external;
    };

} // namespace Hammer

#endif
