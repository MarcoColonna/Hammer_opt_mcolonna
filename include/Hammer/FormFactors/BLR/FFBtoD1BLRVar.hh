///
/// @file  FFBtoD1BLRVar.hh
/// @brief \f$ B \rightarrow D_1 \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD1BLR_VAR
#define HAMMER_FF_BDSSD1BLR_VAR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD1BLRVar final : public FFBLRBase {

    public:

        FFBtoD1BLRVar();

        FFBtoD1BLRVar(const FFBtoD1BLRVar& other) = default;
        FFBtoD1BLRVar& operator=(const FFBtoD1BLRVar& other) = delete;
        FFBtoD1BLRVar(FFBtoD1BLRVar&& other) = delete;
        FFBtoD1BLRVar& operator=(FFBtoD1BLRVar&& other) = delete;
        ~FFBtoD1BLRVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
