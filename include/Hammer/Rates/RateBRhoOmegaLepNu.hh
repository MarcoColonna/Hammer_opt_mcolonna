///
/// @file  RateBRhoOmegaLepNu.hh
/// @brief \f$ B \rightarrow \rho \tau\nu \f$ total rate
/// @brief \f$ B \rightarrow \omega \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_BRHOOMEGALEPNU
#define HAMMER_RATE_BRHOOMEGALEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateBRhoOmegaLepNu final : public RateBase {

    public:

        RateBRhoOmegaLepNu();

        ~RateBRhoOmegaLepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
