///
/// @file  FFBtoDstarBLPRVar.hh
/// @brief \f$ B \rightarrow D^* \f$ BLPRVar form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARBLPRVAR
#define HAMMER_FF_BDSTARBLPRVAR

#include "Hammer/FormFactors/BLPR/FFBLPRVarBase.hh"

namespace Hammer {

    class FFBtoDstarBLPRVar final : public FFBLPRVarBase {

    public:

        FFBtoDstarBLPRVar();

        FFBtoDstarBLPRVar(const FFBtoDstarBLPRVar& other) = default;
        FFBtoDstarBLPRVar& operator=(const FFBtoDstarBLPRVar& other) = delete;
        FFBtoDstarBLPRVar(FFBtoDstarBLPRVar&& other) = delete;
        FFBtoDstarBLPRVar& operator=(FFBtoDstarBLPRVar&& other) = delete;
        ~FFBtoDstarBLPRVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
