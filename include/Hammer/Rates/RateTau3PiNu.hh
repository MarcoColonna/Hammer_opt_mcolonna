///
/// @file  RateTau3PiNu.hh
/// @brief \f$ \tau \rightarrow 3\pi \nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_TAU3PINU
#define HAMMER_RATE_TAU3PINU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateTau3PiNu final : public RateBase {

    public:

        RateTau3PiNu();

        ~RateTau3PiNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
