///
/// @file  AmplBDLepNu.hh
/// @brief \f$ B \rightarrow D \tau\nu \f$ amplitude
/// @brief Also: \f$ B \rightarrow \pi \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_BDLEPNU
#define HAMMER_AMPL_BDLEPNU

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"

namespace Hammer {

    class AmplBDLepNu final : public AmplBToQLepNuBase {

    public:

        AmplBDLepNu();

        ~AmplBDLepNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void addRefs() const final;
    };

} // namespace Hammer

#endif
