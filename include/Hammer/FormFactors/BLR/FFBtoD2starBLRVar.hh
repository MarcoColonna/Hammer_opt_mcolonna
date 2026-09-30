///
/// @file  FFBtoD2starBLRVar.hh
/// @brief \f$ B \rightarrow D_2^* \f$ BLR form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSSD2STARBLR_VAR
#define HAMMER_FF_BDSSD2STARBLR_VAR

#include "Hammer/FormFactors/BLR/FFBLRBase.hh"

namespace Hammer {

    class FFBtoD2starBLRVar final : public FFBLRBase {

    public:

        FFBtoD2starBLRVar();

        FFBtoD2starBLRVar(const FFBtoD2starBLRVar& other) = default;
        FFBtoD2starBLRVar& operator=(const FFBtoD2starBLRVar& other) = delete;
        FFBtoD2starBLRVar(FFBtoD2starBLRVar&& other) = delete;
        FFBtoD2starBLRVar& operator=(FFBtoD2starBLRVar&& other) = delete;
        ~FFBtoD2starBLRVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
