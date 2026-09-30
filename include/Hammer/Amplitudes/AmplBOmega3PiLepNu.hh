///
/// @file  AmplBOmega3PiLepNu.hh
/// @brief \f$ B \rightarrow \omega \tau\nu, \omega \rightarrow 3\pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_BROMEGA3PILEPNU
#define HAMMER_AMPL_BROMEGA3PILEPNU

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"

namespace Hammer {

    class AmplBOmega3PiLepNu final : public AmplBToQLepNuBase {

    public:

        AmplBOmega3PiLepNu();

        ~AmplBOmega3PiLepNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;
    };

} // namespace Hammer

#endif
