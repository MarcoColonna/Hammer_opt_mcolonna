///
/// @file  FFD2starDstarPiPW.hh
/// @brief \f$ D_2^* \rightarrow D^* \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_D2STARDSTARPIPW
#define HAMMER_FF_D2STARDSTARPIPW

#include "Hammer/FormFactors/PW/FFPWBase.hh"

namespace Hammer {

    class FFD2starDstarPiPW final : public FFPWBase {

    public:

        FFD2starDstarPiPW();

        FFD2starDstarPiPW(const FFD2starDstarPiPW& other) = default;
        FFD2starDstarPiPW& operator=(const FFD2starDstarPiPW& other) = delete;
        FFD2starDstarPiPW(FFD2starDstarPiPW&& other) = delete;
        FFD2starDstarPiPW& operator=(FFD2starDstarPiPW&& other) = delete;
        ~FFD2starDstarPiPW() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
