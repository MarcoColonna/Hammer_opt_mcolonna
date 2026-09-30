///
/// @file  FFBtoD1starBLRVar.hh
/// @brief \f$ B \rightarrow D_1^* \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD1STARBLR_VAR
#define HAMMER_FF_BDSSD1STARBLR_VAR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD1starBLRVar final : public FFBLRBase {

    public:

        FFBtoD1starBLRVar();

        FFBtoD1starBLRVar(const FFBtoD1starBLRVar& other) = default;
        FFBtoD1starBLRVar& operator=(const FFBtoD1starBLRVar& other) = delete;
        FFBtoD1starBLRVar(FFBtoD1starBLRVar&& other) = delete;
        FFBtoD1starBLRVar& operator=(FFBtoD1starBLRVar&& other) = delete;
        ~FFBtoD1starBLRVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
