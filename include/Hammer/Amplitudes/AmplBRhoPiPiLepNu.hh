///
/// @file  AmplBRhoPiPiLepNu.hh
/// @brief \f$ B \rightarrow \rho \tau\nu, \rho \rightarrow \pi \pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_BRHOPIPILEPNU
#define HAMMER_AMPL_BRHOPIPILEPNU

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"

namespace Hammer {

    class AmplBRhoPiPiLepNu final : public AmplBToQLepNuBase {

    public:

        AmplBRhoPiPiLepNu();

        ~AmplBRhoPiPiLepNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;
    };

} // namespace Hammer

#endif
