///
/// @file  FFBtoPiBCL.hh
/// @brief \f$ B \rightarrow \pi \f$ BCL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BPIBCL
#define HAMMER_FF_BPIBCL

#include "Hammer/FormFactors/BCL/FFBCLBase.hh"

namespace Hammer {

    class FFBtoPiBCL final : public FFBCLBase {

    public:

        FFBtoPiBCL();

        FFBtoPiBCL(const FFBtoPiBCL& other) = default;
        FFBtoPiBCL& operator=(const FFBtoPiBCL& other) = delete;
        FFBtoPiBCL(FFBtoPiBCL&& other) = delete;
        FFBtoPiBCL& operator=(FFBtoPiBCL&& other) = delete;
        ~FFBtoPiBCL() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
