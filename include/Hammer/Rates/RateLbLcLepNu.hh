///
/// @file  RateLbLcLepNu.hh
/// @brief \f$ \Lambda_b \rightarrow \Lambda_c \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_LBLCLEPNU
#define HAMMER_RATE_LBLCLEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateLbLcLepNu final : public RateBase {

    public:

        RateLbLcLepNu();

        ~RateLbLcLepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
