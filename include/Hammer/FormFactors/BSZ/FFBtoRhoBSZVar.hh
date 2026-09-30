///
/// @file  FFBtoRhoBSZVar.hh
/// @brief \f$ B \rightarrow \rho \f$ BSZ form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BRHOBSZVAR
#define HAMMER_FF_BRHOBSZVAR

#include "Hammer/FormFactors/BSZ/FFBSZBase.hh"

namespace Hammer {

    class FFBtoRhoBSZVar final : public FFBSZBase {

    public:

        FFBtoRhoBSZVar();

        FFBtoRhoBSZVar(const FFBtoRhoBSZVar& other) = default;
        FFBtoRhoBSZVar& operator=(const FFBtoRhoBSZVar& other) = delete;
        FFBtoRhoBSZVar(FFBtoRhoBSZVar&& other) = delete;
        FFBtoRhoBSZVar& operator=(FFBtoRhoBSZVar&& other) = delete;
        ~FFBtoRhoBSZVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
