///
/// @file  RateBD1starLepNu.hh
/// @brief \f$ B \rightarrow D_1^* \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_BDSSD1STARLEPNU
#define HAMMER_RATE_BDSSD1STARLEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateBD1starLepNu final : public RateBase {

    public:

        RateBD1starLepNu();

        ~RateBD1starLepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
