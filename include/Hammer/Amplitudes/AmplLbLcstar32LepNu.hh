///
/// @file  AmplLbLcstar32LepNu.hh
/// @brief \f$ \Lambda_b \rightarrow \Lambda_c^*(2625) \tau\nu \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_LBLCSTAR32LEPNU
#define HAMMER_AMPL_LBLCSTAR32LEPNU

#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"

namespace Hammer {

    class AmplLbLcstar32LepNu final : public AmplBToQLepNuBase {

    public:

        AmplLbLcstar32LepNu();

        ~AmplLbLcstar32LepNu() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void addRefs() const final;
    };

} // namespace Hammer

#endif
