///
/// @file  FFBtoPiBCLVar.hh
/// @brief \f$ B \rightarrow \pi \f$ BCL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BPIBCL_VAR
#define HAMMER_FF_BPIBCL_VAR

#include "Hammer/FormFactors/BCL/FFBCLBase.hh"

namespace Hammer {

    class FFBtoPiBCLVar final : public FFBCLBase {

    public:

        FFBtoPiBCLVar();

        FFBtoPiBCLVar(const FFBtoPiBCLVar& other) = default;
        FFBtoPiBCLVar& operator=(const FFBtoPiBCLVar& other) = delete;
        FFBtoPiBCLVar(FFBtoPiBCLVar&& other) = delete;
        FFBtoPiBCLVar& operator=(FFBtoPiBCLVar&& other) = delete;
        ~FFBtoPiBCLVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
