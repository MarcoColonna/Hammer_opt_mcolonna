///
/// @file  AmplBcJpsiEllEllLepNu.hh
/// @brief \f$ B_c \rightarrow J/\psi \tau\nu, J/\psi \rightarrow \ell \ell \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_BCJPSIELLELLLEPNU
#define HAMMER_AMPL_BCJPSIELLELLLEPNU

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"

namespace Hammer {

    class AmplBcJpsiEllEllLepNu final : public AmplBToQLepNuBase {

    public:

        AmplBcJpsiEllEllLepNu();

        ~AmplBcJpsiEllEllLepNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;
    };

} // namespace Hammer

#endif
