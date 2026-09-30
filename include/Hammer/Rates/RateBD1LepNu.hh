///
/// @file  RateBD1LepNu.hh
/// @brief \f$ B \rightarrow D_1 \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_BDSSD1LEPNU
#define HAMMER_RATE_BDSSD1LEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateBD1LepNu final : public RateBase {

    public:

        RateBD1LepNu();

        ~RateBD1LepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
