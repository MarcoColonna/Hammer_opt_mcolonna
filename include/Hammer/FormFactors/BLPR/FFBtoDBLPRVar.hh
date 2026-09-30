///
/// @file  FFBtoDBLPRVar.hh
/// @brief \f$ B \rightarrow D \f$ BLPRVar form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDBLPRVAR
#define HAMMER_FF_BDBLPRVAR

#include "Hammer/FormFactors/BLPR/FFBLPRVarBase.hh"

namespace Hammer {

    class FFBtoDBLPRVar final : public FFBLPRVarBase {

    public:

        FFBtoDBLPRVar();

        FFBtoDBLPRVar(const FFBtoDBLPRVar& other) = default;
        FFBtoDBLPRVar& operator=(const FFBtoDBLPRVar& other) = delete;
        FFBtoDBLPRVar(FFBtoDBLPRVar&& other) = delete;
        FFBtoDBLPRVar& operator=(FFBtoDBLPRVar&& other) = delete;
        ~FFBtoDBLPRVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
