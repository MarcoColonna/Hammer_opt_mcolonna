///
/// @file  RateBDstarLepNu.hh
/// @brief \f$ B \rightarrow D^* \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_BDSTARLEPNU
#define HAMMER_RATE_BDSTARLEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateBDstarLepNu final : public RateBase {

    public:

        RateBDstarLepNu();

        ~RateBDstarLepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
