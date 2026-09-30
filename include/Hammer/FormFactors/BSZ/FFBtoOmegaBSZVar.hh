///
/// @file  FFBtoOmegaBSZVar.hh
/// @brief \f$ B \rightarrow \omega \f$ BSZ form factors
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BOMEGABSZVAR
#define HAMMER_FF_BOMEGABSZVAR

#include "Hammer/FormFactors/BSZ/FFBSZBase.hh"

namespace Hammer {

    class FFBtoOmegaBSZVar final : public FFBSZBase {

    public:

        FFBtoOmegaBSZVar();

        FFBtoOmegaBSZVar(const FFBtoOmegaBSZVar& other) = default;
        FFBtoOmegaBSZVar& operator=(const FFBtoOmegaBSZVar& other) = delete;
        FFBtoOmegaBSZVar(FFBtoOmegaBSZVar&& other) = delete;
        FFBtoOmegaBSZVar& operator=(FFBtoOmegaBSZVar&& other) = delete;
        ~FFBtoOmegaBSZVar() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
