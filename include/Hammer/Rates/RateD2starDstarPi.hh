///
/// @file  RateD2starDstarPi.hh
/// @brief \f$ D_2^* \rightarrow D^* \pi \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_D2STARDSTARPI
#define HAMMER_RATE_D2STARDSTARPI

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateD2starDstarPi final : public RateBase {

    public:

        RateD2starDstarPi();

        ~RateD2starDstarPi() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
