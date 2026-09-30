///
/// @file  AmplBDstarDPiLepNu.hh
/// @brief \f$ B \rightarrow D^* \tau\nu, D^* \rightarrow D \pi \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_BDSTARDPILEPNU
#define HAMMER_AMPL_BDSTARDPILEPNU

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"

namespace Hammer {

    class AmplBDstarDPiLepNu final : public AmplBToQLepNuBase {

    public:

        AmplBDstarDPiLepNu();

        ~AmplBDstarDPiLepNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void addRefs() const final;
    };

} // namespace Hammer

#endif
