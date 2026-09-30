///
/// @file  AmplTauEllNuNu.hh
/// @brief \f$ \tau-> \ell\nu\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_TAUELLNUNU
#define HAMMER_AMPL_TAUELLNUNU

#include "Hammer/AmplitudeBase.hh"

namespace Hammer {

    class AmplTauEllNuNu final : public AmplitudeBase {

    public:

        AmplTauEllNuNu();

        ~AmplTauEllNuNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void defineSettings() final;
    };

} // namespace Hammer

#endif
