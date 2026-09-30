///
/// @file  FFBstoKBCL.hh
/// @brief \f$ B_s \rightarrow K \f$ BCL form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BSKBCL
#define HAMMER_FF_BSKBCL

#include "Hammer/FormFactors/BCL/FFBCLBase.hh"

namespace Hammer {

    class FFBstoKBCL final : public FFBCLBase {

    public:

        FFBstoKBCL();

        FFBstoKBCL(const FFBstoKBCL& other) = default;
        FFBstoKBCL& operator=(const FFBstoKBCL& other) = delete;
        FFBstoKBCL(FFBstoKBCL&& other) = delete;
        FFBstoKBCL& operator=(FFBstoKBCL&& other) = delete;
        ~FFBstoKBCL() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
