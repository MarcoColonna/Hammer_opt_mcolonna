///
/// @file  AmplLbLcLepNu.hh
/// @brief \f$ \Lambda_b \rightarrow \Lambda_c \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_LBLCLEPNU
#define HAMMER_AMPL_LBLCLEPNU

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"

namespace Hammer {

    class AmplLbLcLepNu final : public AmplBToQLepNuBase {

    public:

        AmplLbLcLepNu();

        ~AmplLbLcLepNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void addRefs() const final;
    };

} // namespace Hammer

#endif
