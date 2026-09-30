///
/// @file  AmplD2starDstarDGamPi.hh
/// @brief \f$ D_2^* \rightarrow D^* \pi, D^* \rightarrow D \gamma \f$ amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_AMPL_D2STARDSTARDGAMMAPI
#define HAMMER_AMPL_D2STARDSTARDGAMMAPI

#include "Hammer/AmplitudeBase.hh"

namespace Hammer {

    class AmplD2starDstarDGamPi final : public AmplitudeBase {

    public:

        AmplD2starDstarDGamPi();

        ~AmplD2starDstarDGamPi() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

    protected:

        void defineSettings() final;
    };

} // namespace Hammer

#endif
