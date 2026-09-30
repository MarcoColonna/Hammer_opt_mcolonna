///
/// @file  FFBtoDBLPRXPVar.hh
/// @brief \f$ B \rightarrow D \f$ BLPRXPVar form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDBLPRXPVAR
#define HAMMER_FF_BDBLPRXPVAR

#include "Hammer/FormFactors/BLPRXP/FFBLPRXPVarBase.hh"

namespace Hammer {

    class FFBtoDBLPRXPVar final : public FFBLPRXPVarBase {

    public:

        FFBtoDBLPRXPVar();

        FFBtoDBLPRXPVar(const FFBtoDBLPRXPVar& other) = default;
        FFBtoDBLPRXPVar& operator=(const FFBtoDBLPRXPVar& other) = delete;
        FFBtoDBLPRXPVar(FFBtoDBLPRXPVar&& other) = delete;
        FFBtoDBLPRXPVar& operator=(FFBtoDBLPRXPVar&& other) = delete;
        ~FFBtoDBLPRXPVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
