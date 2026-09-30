///
/// @file  AmplD1DstarDGamPi.hh
/// @brief \f$ D_1 \rightarrow D^* \pi, D^* \rightarrow D \gamma \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_D1DSTARDGAMMAPI
#define HAMMER_AMPL_D1DSTARDGAMMAPI

#include "Hammer/AmplitudeBase.hh"

namespace Hammer {

    class AmplD1DstarDGamPi final : public AmplitudeBase {

    public:

        AmplD1DstarDGamPi();

        ~AmplD1DstarDGamPi() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void defineSettings() final;
    };

} // namespace Hammer

#endif
