///
/// @file  FFBstoKBCLVar.hh
/// @brief \f$ B_s \rightarrow K \f$ BCL form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BSKBCL_VAR
#define HAMMER_FF_BSKBCL_VAR

#include "Hammer/FormFactors/BCL/FFBCLBase.hh"

namespace Hammer {

    class FFBstoKBCLVar final : public FFBCLBase {

    public:

        FFBstoKBCLVar();

        FFBstoKBCLVar(const FFBstoKBCLVar& other) = default;
        FFBstoKBCLVar& operator=(const FFBstoKBCLVar& other) = delete;
        FFBstoKBCLVar(FFBstoKBCLVar&& other) = delete;
        FFBstoKBCLVar& operator=(FFBstoKBCLVar&& other) = delete;
        ~FFBstoKBCLVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) override;

        void defineSettings() final;

        void addRefs() const final;
    };

} // namespace Hammer

#endif
