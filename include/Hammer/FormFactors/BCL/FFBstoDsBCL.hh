///
/// @file  FFBstoDsBCL.hh
/// @brief \f$ B_s \rightarrow D_s \f$ BCL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BSDSBCL
#define HAMMER_FF_BSDSBCL

#include "Hammer/FormFactors/BCL/FFBCLBase.hh"

namespace Hammer {

    class FFBstoDsBCL final : public FFBCLBase {

    public:

        FFBstoDsBCL();

        FFBstoDsBCL(const FFBstoDsBCL& other) = default;
        FFBstoDsBCL& operator=(const FFBstoDsBCL& other) = delete;
        FFBstoDsBCL(FFBstoDsBCL&& other) = delete;
        FFBstoDsBCL& operator=(FFBstoDsBCL&& other) = delete;
        ~FFBstoDsBCL() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
