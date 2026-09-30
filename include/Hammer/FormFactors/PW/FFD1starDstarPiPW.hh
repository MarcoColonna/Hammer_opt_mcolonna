///
/// @file  FFD1starDstarPiPW.hh
/// @brief \f$ D_1^* \rightarrow D^* \pi \f$ partial wave coefficients
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FF_D1STARDSTARPIPW
#define HAMMER_FF_D1STARDSTARPIPW

#include "Hammer/FormFactors/PW/FFPWBase.hh"

namespace Hammer {

    class FFD1starDstarPiPW final : public FFPWBase {

    public:

        FFD1starDstarPiPW();

        FFD1starDstarPiPW(const FFD1starDstarPiPW& other) = default;
        FFD1starDstarPiPW& operator=(const FFD1starDstarPiPW& other) = delete;
        FFD1starDstarPiPW(FFD1starDstarPiPW&& other) = delete;
        FFD1starDstarPiPW& operator=(FFD1starDstarPiPW&& other) = delete;
        ~FFD1starDstarPiPW() final = default;

        void eval(const Particle& parent, const ParticleList& daughters, const ParticleList& references) final;

        std::unique_ptr<FormFactorBase> clone(const std::string& label) final;

    protected:

        void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses = {}) final;

        void defineSettings() final;
    };

} // namespace Hammer

#endif
