///
/// @file  AmplD1DstarDPiPi.hh
/// @brief \f$ D_1 \rightarrow D^* \pi, D^* \rightarrow D \pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_D1DSTARDPIPI
#define HAMMER_AMPL_D1DSTARDPIPI

#include "Hammer/AmplitudeBase.hh"

namespace Hammer {

    class AmplD1DstarDPiPi final : public AmplitudeBase {

    public:

        AmplD1DstarDPiPi();

        ~AmplD1DstarDPiPi() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void defineSettings() final;
    };

} // namespace Hammer

#endif
