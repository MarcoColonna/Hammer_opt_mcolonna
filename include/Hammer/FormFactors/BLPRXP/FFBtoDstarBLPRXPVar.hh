///
/// @file  FFBtoDstarBLPRXPVar.hh
/// @brief \f$ B \rightarrow D^* \f$ BLPRXPVar form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBLPRXPVAR
#define HAMMER_FF_BDSTARBLPRXPVAR

#include "Hammer/FormFactors/BLPRXP/FFBLPRXPVarBase.hh"

namespace Hammer {

    class FFBtoDstarBLPRXPVar final : public FFBLPRXPVarBase {

    public:

        FFBtoDstarBLPRXPVar();

        FFBtoDstarBLPRXPVar(const FFBtoDstarBLPRXPVar& other) = default;
        FFBtoDstarBLPRXPVar& operator=(const FFBtoDstarBLPRXPVar& other) = delete;
        FFBtoDstarBLPRXPVar(FFBtoDstarBLPRXPVar&& other) = delete;
        FFBtoDstarBLPRXPVar& operator=(FFBtoDstarBLPRXPVar&& other) = delete;
        ~FFBtoDstarBLPRXPVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
