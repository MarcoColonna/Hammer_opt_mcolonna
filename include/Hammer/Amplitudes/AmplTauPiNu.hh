///
/// @file  AmplTauPiNu.hh
/// @brief \f$ \tau-> \pi\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_TAUPINU
#define HAMMER_AMPL_TAUPINU

#include "Hammer/AmplitudeBase.hh"

namespace Hammer {

    class AmplTauPiNu final : public AmplitudeBase {

    public:

        AmplTauPiNu();

        ~AmplTauPiNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void defineSettings() final;
    };

} // namespace Hammer

#endif
