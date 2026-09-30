///
/// @file  RateBDLepNu.hh
/// @brief \f$ B \rightarrow D \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_BDLEPNU
#define HAMMER_RATE_BDLEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateBDLepNu final : public RateBase {

    public:

        RateBDLepNu();

        ~RateBDLepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
