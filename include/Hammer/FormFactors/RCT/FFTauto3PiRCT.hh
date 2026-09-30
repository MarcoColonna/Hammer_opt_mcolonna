///
/// @file  FFTauto3PiRCT.hh
/// @brief \f$ \tau^+ \rightarrow \bar\nu\pi^+\pi^+\pi^- \f$ form factors
/// see 1203.3955 and 1310.1053
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_TAU3PIRCT
#define HAMMER_FF_TAU3PIRCT

#include "Hammer/FormFactors/RCT/FFRCTBase.hh"

namespace Hammer {

    class FFTauto3PiRCT final : public FFRCTBase {

    public:

        FFTauto3PiRCT();

        FFTauto3PiRCT(const FFTauto3PiRCT& other) = default;
        FFTauto3PiRCT& operator=(const FFTauto3PiRCT& other) = delete;
        FFTauto3PiRCT(FFTauto3PiRCT&& other) = delete;
        FFTauto3PiRCT& operator=(FFTauto3PiRCT&& other) = delete;
        ~FFTauto3PiRCT() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
