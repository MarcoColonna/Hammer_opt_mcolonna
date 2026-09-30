///
/// @file  FFBtoDCLNVar.hh
/// @brief \f$ B \rightarrow D \f$ CLN form factors
/// @brief Ported from EvtGen HQET3 (custom class; F Bernlochner)
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDCLNVAR
#define HAMMER_FF_BDCLNVAR

#include "Hammer/FormFactors/CLN/FFCLNBase.hh"

namespace Hammer {

    class FFBtoDCLNVar final : public FFCLNBase {

    public:

        FFBtoDCLNVar();

        FFBtoDCLNVar(const FFBtoDCLNVar& other) = default;
        FFBtoDCLNVar& operator=(const FFBtoDCLNVar& other) = delete;
        FFBtoDCLNVar(FFBtoDCLNVar&& other) = delete;
        FFBtoDCLNVar& operator=(FFBtoDCLNVar&& other) = delete;
        ~FFBtoDCLNVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
