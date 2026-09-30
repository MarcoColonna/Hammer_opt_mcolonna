///
/// @file  RateD2starDPi.hh
/// @brief \f$ D_2^* \rightarrow D \pi \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_D2STARDPI
#define HAMMER_RATE_D2STARDPI

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateD2starDPi final : public RateBase {

    public:

        RateD2starDPi();

        ~RateD2starDPi() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
