///
/// @file  FFBtoDstarCLNVar.hh
/// @brief \f$ B \rightarrow D^* \f$ CLN form factors with variations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BDSTARCLNVAR
#define HAMMER_FF_BDSTARCLNVAR

#include "Hammer/FormFactors/CLN/FFCLNBase.hh"

namespace Hammer {

    class FFBtoDstarCLNVar final : public FFCLNBase {

    public:

        FFBtoDstarCLNVar();

        FFBtoDstarCLNVar(const FFBtoDstarCLNVar& other) = default;
        FFBtoDstarCLNVar& operator=(const FFBtoDstarCLNVar& other) = delete;
        FFBtoDstarCLNVar(FFBtoDstarCLNVar&& other) = delete;
        FFBtoDstarCLNVar& operator=(FFBtoDstarCLNVar&& other) = delete;
        ~FFBtoDstarCLNVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
