///
/// @file  FFBtoPiGKvD.hh
/// @brief \f$ B \rightarrow \pi \f$ GKvD form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BPIGKVD
#define HAMMER_FF_BPIGKVD

#include "Hammer/FormFactors/BSZ/FFBSZBase.hh"

namespace Hammer {

    class FFBtoPiGKvD final : public FFBSZBase {

    public:

        FFBtoPiGKvD();

        FFBtoPiGKvD(const FFBtoPiGKvD& other) = default;
        FFBtoPiGKvD& operator=(const FFBtoPiGKvD& other) = delete;
        FFBtoPiGKvD(FFBtoPiGKvD&& other) = delete;
        FFBtoPiGKvD& operator=(FFBtoPiGKvD&& other) = delete;
        ~FFBtoPiGKvD() final = default;


        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
