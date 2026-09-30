///
/// @file  RateLbLcstar32LepNu.hh
/// @brief \f$ \Lambda_b \rightarrow \Lambda_c^*(2625) \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_LBLCSTAR32LEPNU
#define HAMMER_RATE_LBLCSTAR32LEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateLbLcstar32LepNu final : public RateBase {

    public:

        RateLbLcstar32LepNu();

        ~RateLbLcstar32LepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
