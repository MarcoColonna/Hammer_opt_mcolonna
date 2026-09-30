///
/// @file  FFD2starDPiPW.hh
/// @brief \f$ D_2^* \rightarrow D \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_D2STARDPIPW
#define HAMMER_FF_D2STARDPIPW

#include "Hammer/FormFactors/PW/FFPWBase.hh"

namespace Hammer {

    class FFD2starDPiPW final : public FFPWBase {

    public:

        FFD2starDPiPW();

        FFD2starDPiPW(const FFD2starDPiPW& other) = default;
        FFD2starDPiPW& operator=(const FFD2starDPiPW& other) = delete;
        FFD2starDPiPW(FFD2starDPiPW&& other) = delete;
        FFD2starDPiPW& operator=(FFD2starDPiPW&& other) = delete;
        ~FFD2starDPiPW() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
