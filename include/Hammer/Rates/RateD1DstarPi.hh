///
/// @file  RateD1DstarPi.hh
/// @brief \f$ D_1 \rightarrow D^* \pi \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_D1DSTARPI
#define HAMMER_RATE_D1DSTARPI

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateD1DstarPi final : public RateBase {

    public:

        RateD1DstarPi();

        ~RateD1DstarPi() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
