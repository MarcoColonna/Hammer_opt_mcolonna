///
/// @file  FFBtoD0starBLRVar.hh
/// @brief \f$ B \rightarrow D_0^* \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD0STARBLR_VAR
#define HAMMER_FF_BDSSD0STARBLR_VAR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD0starBLRVar final : public FFBLRBase {

    public:

        FFBtoD0starBLRVar();

        FFBtoD0starBLRVar(const FFBtoD0starBLRVar& other) = default;
        FFBtoD0starBLRVar& operator=(const FFBtoD0starBLRVar& other) = delete;
        FFBtoD0starBLRVar(FFBtoD0starBLRVar&& other) = delete;
        FFBtoD0starBLRVar& operator=(FFBtoD0starBLRVar&& other) = delete;
        ~FFBtoD0starBLRVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
