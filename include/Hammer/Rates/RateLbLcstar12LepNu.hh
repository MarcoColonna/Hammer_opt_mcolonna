///
/// @file  RateLbLcstar12LepNu.hh
/// @brief \f$ \Lambda_b \rightarrow \Lambda_c^*(2595) \tau\nu \f$ total rate
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_RATE_LBLCSTAR12LEPNU
#define HAMMER_RATE_LBLCSTAR12LEPNU

#include "Hammer/RateBase.hh"

namespace Hammer {

    class RateLbLcstar12LepNu final : public RateBase {

    public:

        RateLbLcstar12LepNu();

        ~RateLbLcstar12LepNu() final = default;

    protected:

        Tensor evalAtPSPoint(const std::vector<double>& point) final;
    };

} // namespace Hammer

#endif
