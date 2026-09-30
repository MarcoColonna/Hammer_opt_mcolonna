///
/// @file  FFBtoRhoOmegaISGW2.hh
/// @brief \f$ B \rightarrow \rho/\omega \f$ ISGW2 form factors
/// @brief Ported directly from EvtGen
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_BRHOOMEGAISGW2
#define HAMMER_FF_BRHOOMEGAISGW2

#include "Hammer/FormFactors/ISGW2/FFISGW2Base.hh"

namespace Hammer {

    class FFBtoRhoOmegaISGW2 final : public FFISGW2Base {

    public:

        FFBtoRhoOmegaISGW2();

        FFBtoRhoOmegaISGW2(const FFBtoRhoOmegaISGW2& other) = default;
        FFBtoRhoOmegaISGW2& operator=(const FFBtoRhoOmegaISGW2& other) = delete;
        FFBtoRhoOmegaISGW2(FFBtoRhoOmegaISGW2&& other) = delete;
        FFBtoRhoOmegaISGW2& operator=(FFBtoRhoOmegaISGW2&& other) = delete;
        ~FFBtoRhoOmegaISGW2() final = default;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
